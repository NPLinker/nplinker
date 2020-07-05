import sys
import logging
import copy

import numpy as np

from .metabolomics import Spectrum
from .metabolomics import MolecularFamily

from .genomics import GCF, BGC

from .config import Config, Args
from .loader import DatasetLoader

from .pickler import save_pickled_data

from .scoring.methods import MetcalfScoring, RosettaScoring, TestScoring
from .scoring.methods import LinkCollection

from .logconfig import LogConfig
logger = LogConfig.getLogger(__file__)


class NPLinker(object):

    # allowable types for objects to be passed to scoring methods
    OBJ_CLASSES = [Spectrum, MolecularFamily, GCF, BGC]
    # default set of enabled scoring methods
    SCORING_METHODS =  {
        MetcalfScoring.NAME: MetcalfScoring,
        TestScoring.NAME: TestScoring,
        RosettaScoring.NAME: RosettaScoring, 
                       }

    def __init__(self, userconfig=None, platform_id=None):
        # TODO update docstring
        """Initialise an NPLinker instance, automatically loading a dataset and generating scores.

        NPLinker instances can be configured in multiple ways, in ascending order of priority:
            1. A global user-level default configuration file in TOML format, found in the directory:
                    $XDG_CONFIG_HOME/nplinker/nplinker.toml
            2. A local TOML configuration file
            3. Command-line arguments / supplying a manually constructed dict instance

        The user-level configuration file will be created when you first create an NPLinker instance. 
        It contains sensible default values for each setting and is intended to be copied and edited to
        produce dataset-specific configuration files, which will then override any parameters shared
        with the user-level file. To load such a file, simply set the "userconfig" parameter to a 
        string containing the filename. 

        It's also possible to selectively override configuration file parameters by 
        supplying command-line arguments (if running nplinker.py as a script), or by passing
        a dict with a structure corresponding to the configuration file format to this method.

        Some examples may make the various possible combinations a bit clearer:
            # load the default/user-level configuration file and nothing else
            npl = NPLinker()

            # override the default file with a different one
            npl = NPLinker('myconfig.toml')

            # the same thing but running NPLinker as a script
            > python nplinker.py --config "myconfig.toml"

            # use the defaults from the user-level config while modifying the root path
            # to load the dataset from (this is the minimum you would need to change in the
            # default config file)
            npl = NPLinker({'dataset': {'root': '/path/to/dataset'}})

            # the same thing running NPLinker as a script
            > python nplinker.py --dataset.root /path/to/dataset
    
        Args:
            userconfig: supplies user-defined configuration data. May take one of 3 types:
                - None: just load the user-level default configuration file
                - str: treat as filename of a local configuration file to load 
                        (overriding the defaults)
                - dict: contents will be used to override values in the dict generated 
                        from loading the configuration file(s)
        """

        # if userconfig is None => create Config() from empty dict
        # if userconfig is a string => create a dict with 'config' key and string as filename
        # if userconfig is a dict => pass it to Config() directly
        if userconfig is None:
            userconfig = {}
        elif isinstance(userconfig, str):
            userconfig = {'config': userconfig}
        elif not isinstance(userconfig, dict):
            raise Exception('Invalid type for userconfig (should be None/str/dict, found "{}")'.format(type(userconfig)))

        if platform_id is not None:
            logger.debug('Setting project ID to {}'.format(platform_id))
            userconfig['dataset'] = {'platform_id': platform_id }

        self._config = Config(userconfig)

        # configure logging based on the supplied config params
        LogConfig.setLogLevelStr(self._config.config['loglevel'])
        logfile = self._config.config['logfile']
        if len(logfile) > 0:
            logfile_dest = logging.FileHandler(logfile)
            # if we want to log to stdout plus logfile, add the new destination
            if self._config.config.get('log_to_stdout', True): # default to True
                LogConfig.addLogDestination(logfile_dest)
            else:
                # otherwise overwrite the default stdout destination
                LogConfig.setLogDestination(logfile_dest)

        # this object takes care of figuring out the locations of all the relevant files/folders
        # and will show error/warning if any missing (depends if optional or not)
        self._loader = DatasetLoader(self._config.config)
    
        self._spectra = []
        self._bgcs = []
        self._gcfs = []
        self._strains = None
        self._metadata = {}
        self._molfams = []
        self._mibig_bgc_dict = {}

        self._bgc_lookup = {}
        self._spec_lookup = {}

        self._scoring_methods = {}
        config_methods = self._config.config.get('scoring_methods', [])
        for name, method in NPLinker.SCORING_METHODS.items():
            if len(config_methods) == 0 or name in config_methods:
                self._scoring_methods[name] = method
                logger.debug('Enabled scoring method: {}'.format(name))

        self._scoring_methods_setup_complete = {name: False for name in self._scoring_methods.keys()}

        self._datalinks = None 

        self._repro_data = {}
        repro_file = self._config.config['repro_file']
        if len(repro_file) > 0:
            self._save_repro_data(repro_file)

    def _collect_repro_data(self):
        """Creates a dict containing data to aid reproducible runs of nplinker.

        This method creates a dict containing various bits of information about
        the current execution of nplinker. This data will typically be saved to
        a file in order to aid reproducibility using :func:`_save_repro_data`. 

        TODO describe contents

        Returns:
            A dict containing the information described above
        """

        self._repro_data = {}
        # TODO best way to embed latest git commit hash? probably with a packaging script...
        # TODO versions of all Python dependencies used (just rely on Pipfile.lock here?)
        
        # insert command line arguments
        self._repro_data['args'] = {}
        for i, arg in enumerate(sys.argv):
            self._repro_data['args'][i] = arg

        # TODO anything else to include here?

        return self._repro_data

    def _save_repro_data(self, filename):
        self._collect_repro_data()
        with open(filename, 'wb') as repro_file:
            # TODO is pickle the best format to use?
            save_pickled_data(self._repro_data, repro_file)
            logger.info('Saving reproducibility data to {}'.format(filename))

    @property
    def config(self):
        """Returns a copy of the data parsed from the configuration file as a dict"""
        return copy.deepcopy(self._config.config)

    @property
    def root_dir(self):
        """Returns path to nplinker external dataset directory (user-configured)"""
        return self._loader._root

    @property
    def dataset_id(self):
        """Returns dataset "ID". For local datasets this will just be the last component
        of the directory path, but for platform datasets it will be the platform 
        project ID"""
        return self._loader.dataset_id

    @property
    def data_dir(self):
        """Returns path to nplinker/data directory (files packaged with the app itself)"""
        return self._loader.datadir

    @property
    def gnps_params(self):
        """Returns a dict containing data from GNPS params.xml (if available)"""
        return self._loader.gnps_params

    @property
    def dataset_description(self):
        """Returns a string containing the content of any supplied description.txt file"""
        return self._loader.description_text

    @property
    def bigscape_cutoff(self):
        """Returns the current BiGSCAPE clustering cutoff value"""
        return self._loader._bigscape_cutoff

    def load_data(self, new_bigscape_cutoff=None, met_only=False):
        """Loads the basic components of a dataset.

        This method is responsible for loading the various pieces of the supplied dataset into
        memory and doing any initial parsing/object manipulation required. After it completes,
        applications can access the lists of GCFs, Spectra, MolecularFamilies and strains 
        using the corresponding properties of the NPLinker class.

        Returns:
            True if successful, False otherwise
        """
        # typical case where load_data is being called with no params
        if new_bigscape_cutoff is None:
            logger.debug('load_data (normal case, full load, met_only={})'.format(met_only))
            self._loader.validate()

            if not self._loader.load(met_only=met_only):
                return False
        else:
            logger.debug('load_data with new cutoff = {}'.format(new_bigscape_cutoff))
            # 1. change the cutoff (which by itself doesn't do anything)
            self._loader._bigscape_cutoff = new_bigscape_cutoff
            # 2. reload the strain mappings (MiBIG filtering may have removed strains
            # that were originally present, need to restore them all so the filtering
            # won't break when it runs again in next stage)
            self._loader._load_strain_mappings()
            # 3. reload the genomics data with the new cutoff applied
            self._loader._load_genomics()

        self._spectra = self._loader.spectra
        self._molfams = self._loader.molfams
        self._bgcs = self._loader.bgcs
        self._gcfs = self._loader.gcfs
        self._mibig_bgc_dict = self._loader.mibig_bgc_dict
        self._strains = self._loader.strains
        self._product_types = self._loader.product_types

        logger.debug('Generating lookup tables: genomics')
        self._bgc_lookup = {}
        for i, bgc in enumerate(self._bgcs):
            self._bgc_lookup[bgc.name] = i
        
        self._gcf_lookup = {}
        for i, gcf in enumerate(self._gcfs):
            self._gcf_lookup[gcf.gcf_id] = i

        # don't need to do these two if cutoff changed (indicating genomics data
        # was reloaded but not metabolomics)
        if new_bigscape_cutoff is None:
            logger.debug('Generating lookup tables: metabolomics')
            self._spec_lookup = {}
            for i, spec in enumerate(self._spectra):
                self._spec_lookup[spec.spectrum_id] = i

            self._molfam_lookup = {}
            for i, molfam in enumerate(self._molfams):
                self._molfam_lookup[molfam.id] = i

        logger.debug('load_data: completed')
        return True

    def get_links(self, input_objects, scoring_methods, and_mode=True):
        """Find links for a set of input objects (BGCs/GCFs/Spectra/MolFams)
        
        Given a set of input objects and 1 or more NPLinker scoring methods,
        this method will return a dict mapping input objects with links to
        the result(s) of the scoring method(s) selected. 

        The input objects can be any mix of the following types:
            - BGC
            - GCF
            - Spectrum
            - MolecularFamily

        # TODO update once finished below
        
        Returns:
            A dict keyed by objects from the input_objects list, but only containing
            those objects that had links returned. The values of the dict will be
            a set of the results from each of the enabled scoring methods for the
            given object. The content of these results depends on the scoring method. 
        """

        if isinstance(input_objects, list) and len(input_objects) == 0:
            raise Exception('input_objects length must be > 0')

        elif isinstance(scoring_methods, list) and  len(scoring_methods) == 0:
            raise Exception('scoring_methods length must be > 0')

        # check if input_objects is a list of lists. if so there should be one
        # entry for each supplied method for it to be a valid parameter
        if isinstance(input_objects[0], list):
            if len(input_objects) != len(scoring_methods):
                raise Exception('Number of input_objects lists must match number of scoring_methods (found: {}, expected: {})'.format(len(input_objects), len(scoring_methods)))

        # TODO check scoring_methods only contains ScoringMethod-derived instances

        # want everything to be in lists of lists
        if not isinstance(input_objects, list) or (isinstance(input_objects, list) and not isinstance(input_objects[0], list)):
            input_objects = [input_objects]
        if not isinstance(scoring_methods, list):
            scoring_methods = [scoring_methods]

        logger.debug('get_links: {} object sets, {} methods'.format(len(input_objects), len(scoring_methods)))

        # copy the object set if required to make up the numbers
        if len(input_objects) != len(scoring_methods):
            if len(scoring_methods) < len(input_objects):
                raise Exception('Number of scoring methods must be >= number of input object sets')
            elif (len(scoring_methods) > len(input_objects)) and len(input_objects) != 1:
                raise Exception('Mismatch between number of scoring methods and input objects ({} vs {})'.format(len(scoring_methods), len(input_objects)))
            elif len(scoring_methods) > len(input_objects):
                # this is a special case for convenience: pass in 1 set of objects and multiple methods,
                # result is that set is used for all methods
                logger.debug('Duplicating input object set')
                while len(input_objects) < len(scoring_methods):
                    input_objects.append(input_objects[0])
                    logger.debug('Duplicating input object set')

        object_counts = {}
        link_collection = LinkCollection(and_mode)

        for i, method in enumerate(scoring_methods):
            # do any one-off initialisation required by this method
            if not self._scoring_methods_setup_complete[method.name]:
                logger.debug('Doing one-time setup for {}'.format(method.name))
                self._scoring_methods[method.name].setup(self)
                self._scoring_methods_setup_complete[method.name] = True

            # should construct a dict of {object_with_link: <link_data>} entries
            objects_for_method = input_objects[i]
            logger.debug('Calling scoring method {} on {} objects'.format(method.name, len(objects_for_method)))
            link_collection = method.get_links(objects_for_method, link_collection)

        if not self._datalinks:
            logger.debug('Creating internal datalinks object')
            self._datalinks = self.scoring_method(MetcalfScoring.NAME).datalinks
            logger.debug('Created internal datalinks object')

        if len(link_collection) == 0:
            logger.debug('No links found or remaining after merging all method results!')

        # populate shared strain info
        logger.debug('Calculating shared strain information...')
        # TODO more efficient version?
        for source, link_data in link_collection.links.items():
            if isinstance(source, BGC):
                logger.debug('Cannot determine shared strains for BGC input!')
                break

            targets = list(filter(lambda x: not isinstance(x, BGC), link_data.keys()))
            if len(targets) > 0:
                shared_strains = self._datalinks.common_strains([source], targets, True)
                for objpair in shared_strains.keys():
                    shared_strains[objpair] = [self._strains.lookup_index(x) for x in shared_strains[objpair]]

                if isinstance(source, GCF):
                    for target, link in link_data.items():
                        if (target, source) in shared_strains:
                            link.shared_strains = shared_strains[(target, source)]
                else:
                    for target, link in link_data.items():
                        if (source, target) in shared_strains:
                            link.shared_strains = shared_strains[(source, target)]

        logger.debug('Finished calculating shared strain information')

        logger.debug('Final size of link collection is {}'.format(len(link_collection)))
        return link_collection

    def get_common_strains(self, objects_a, objects_b, filter_no_shared=True):
        """
        Allows for retrieval of shared strains for arbitrary pairs of objects.

        Takes two lists of objects as input. Typically one list will be MolecularFamily
        or Spectrum objects and the other GCF (which list is which doesn't matter). The
        return value is a dict mapping pairs of objects to lists of Strain objects 
        shared by that pair. This list may be empty if filter_no_shared is False. If
        filter_no_shared is True, every entry in the dict with no shared strains will
        be removed before it is returned. 

        Args:
            objects_a: a list of Spectrum/MolecularFamily/GCF objects
            objects_b: a list of Spectrum/MolecularFamily/GCF objects 
            filter_no_shared: remove result entries if no shared strains exist

        Returns:
            A dict mapping pairs of objects (obj1, obj2) to lists of Strain objects.

            NOTE: The ordering of the pairs is *fixed* to be (metabolomic, genomic). In
            other words, if objects_a = [GCF1, GC2, ...] and objects_b = [Spectrum1, 
            Spectrum2, ...], the object pairs will be (Spectrum1, GCF1), (Spectrum2,
            GCF2), and so on. The same applies if objects_a and objects_b are swapped,
            the metabolomic objects (Spectrum or MolecularFamily) will be the obj1
            entry in each pair. 
        """
        if not self._datalinks:
            self._datalinks = self.scoring_method(MetcalfScoring.NAME).datalinks

        # this is a dict with structure:
        #   (Spectrum/MolecularFamily, GCF) => list of strain indices
        common_strains = self._datalinks.common_strains(objects_a, objects_b, filter_no_shared)

        # replace the lists of strain indices with actual strain objects
        for objpair in common_strains.keys():
            common_strains[objpair] = [self._strains.lookup_index(x) for x in common_strains[objpair]]
    
        return common_strains

    def has_bgc(self, name):
        return name in self._bgc_lookup

    def lookup_bgc(self, name):
        if name not in self._bgc_lookup:
            return None

        return self._bgcs[self._bgc_lookup[name]]

    def lookup_gcf(self, gcf_id):
        if gcf_id not in self._gcf_lookup:
            return None

        return self._gcfs[self._gcf_lookup[gcf_id]]

    def lookup_spectrum(self, name):
        if name not in self._spec_lookup:
            return None

        return self._spectra[self._spec_lookup[name]]

    @property
    def strains(self):
        return self._strains

    @property
    def bgcs(self):
        return self._bgcs

    @property
    def gcfs(self):
        return self._gcfs

    @property
    def spectra(self):
        return self._spectra

    @property
    def molfams(self):
        return self._molfams

    @property
    def metadata(self):
        return self._metadata

    @property
    def mibig_bgc_dict(self):
        return self._mibig_bgc_dict

    @property
    def product_types(self):
        """
        Returns a list of the available BiGSCAPE product types in current dataset
        """
        return self._product_types

    @property
    def repro_data(self):
        """
        Returns the dict of reproducibility data
        """
        return self._repro_data

    @property
    def scoring_methods(self):
        """
        Returns the set of available scoring method names
        """
        return list(self._scoring_methods.keys())

    def scoring_method(self, name):
        """
        Return an instance of the given scoring method
        """
        if not self._scoring_methods_setup_complete[name]:
            self._scoring_methods[name].setup(self)
            self._scoring_methods_setup_complete[name] = True

        return self._scoring_methods.get(name, None)(self)

if __name__ == "__main__":
    # can set default logging configuration this way...
    LogConfig.setLogLevel(logging.DEBUG)

    # initialise NPLinker from the command-line arguments
    npl = NPLinker(Args().get_args())

    # load the dataset
    if not npl.load_data():
        print('Failed to load the dataset!')
        sys.exit(-1)

    # create a metcalf scoring object
    mc = npl.scoring_method('metcalf')
    # set a scoring cutoff threshold
    mc.cutoff = 3.0

    # pick some GCFs to get links for
    test_gcfs = npl.gcfs[:10]

    # tell nplinker to find links for this set of GCFs using metcalf scoring
    results = npl.get_links(test_gcfs, mc)

    # check if any links were found
    if len(results) == 0:
        print('No links found!')
        sys.exit(0)

    # the "result" object will be a LinkCollection, holding all the information
    # returned by the scoring method(s) used
    print('{} total links found'.format(len(results)))

    # display some information about each object and its links
    for obj, result in results.links.items():
        print('Results for object: {}, {} total links, {} methods used'.format(obj, len(result), results.method_count))

        # get links for this object, sorted by metcalf score
        sorted_links = results.get_sorted_links(mc, obj)
        for link_data in sorted_links:
            print('  --> [{}] {} | {} | # shared_strains = {}'.format(','.join(method.name for method in link_data.methods), link_data.target, link_data[mc], len(link_data.shared_strains)))

