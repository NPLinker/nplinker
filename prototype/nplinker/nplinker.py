import sys
import logging
import pickle

import numpy as np

from .metabolomics import Spectrum
from .metabolomics import MolecularFamily

from .genomics import GCF

from .data_linking import DataLinks, LinkLikelihood
from .data_linking import LinkFinder, SCORING_METHODS

from .config import Config, Args
from .loader import DatasetLoader

from .logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class NPLinker(object):

    OBJ_CLASSES = [Spectrum, MolecularFamily, GCF]
    # enumeration for accessing results of LinkFinder.get_links, which are (3, num_links) arrays:
    # - R_SRC_ID: the ID of an object that was supplied as input to get_links
    # - R_DST_ID: the ID of an object that was discovered to have a link to an input object
    # - R_SCORE: the score for the link between a pair of objects
    R_SRC_ID, R_DST_ID, R_SCORE = range(3)

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
            LogConfig.setLogDestination(logfile_dest)

        # this object takes care of figuring out the locations of all the relevant files/folders
        # and will show error/warning if any missing (depends if optional or not)
        self._loader = DatasetLoader(self._config.config)
    
        self._spectra = []
        self._bgcs = []
        self._gcfs = []
        self._strains = None
        self._metadata = {}
        self._families = []
        self._mibig_bgc_dict = {}

        self._bgc_lookup = {}
        self._spec_lookup = {}

        self._datalinks = None
        self._linkfinder = None

        self._links = {}

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
            pickle.dump(self._repro_data, repro_file)
            logger.info('Saving reproducibility data to {}'.format(filename))

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

        # clear any existing link/score data
        self.clear_links()

        # typical case where load_data is being called with no params
        if new_bigscape_cutoff is None:
            logger.debug('load_data (normal case, full load)')
            self._loader.validate()

            if not self._loader.load(met_only):
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


    def _generate_scores(self, dl, lf):
        # if metcalf scoring enabled
        if self.scoring.metcalf.enabled:
            lf.metcalf_scoring(dl, type='spec-gcf')
            lf.metcalf_scoring(dl, type='fam-gcf')

        # if hg scoring enabled
        if self.scoring.hg.enabled:
            lf.hg_scoring(dl, type='spec-gcf')
            lf.hg_scoring(dl, type='fam-gcf')

        # if likescore scoring enabled
        if self.scoring.likescore.enabled:
            ll = LinkLikelihood()
            ll.calculate_likelihoods(dl, type='spec-gcf')
            ll.calculate_likelihoods(dl, type='fam-gcf')
            lf.likelihood_scoring(dl, ll, type='spec-gcf')
            lf.likelihood_scoring(dl, ll, type='fam-gcf')

    def process_dataset(self):
        """Construct the DataLinks and LinkFinder objects from the loaded dataset.

        Deals with initialising all the objects used for scoring/linking.
        """

        if len(self._spectra) == 0 or len(self._gcfs) == 0 or len(self._strains) == 0:
            logger.info('process_dataset: calling load_data')
            self.load_data()

        self._datalinks = DataLinks()
        self._datalinks.load_data(self._spectra, self._gcfs, self._strains)
        self._datalinks.find_correlations()
        self._linkfinder = LinkFinder()

        # generate the actual scores from the standard DataLinks object
        logger.debug('Generating scores, enabled methods={}'.format(self.scoring.enabled()))
        self._generate_scores(self._datalinks, self._linkfinder)

        logger.debug('Scoring matrices created')

        self.clear_links()
        return True

    def clear_links(self):
        """Clear any previously retrieved links"""
        self._links = {}

    def get_links(self, objects, scoring_method, datalinks=None):
        """Collects sets of links between a given set of objects and their counterparts.

        This method allows an application to pass in one or more objects of type GCF, Spectrum,
        or MolecularFamily, and retrieve the set of links to their counterpart objects (e.g.
        GCF for Spectrum input) that meet the scoring criteria for the selected scoring
        method (see ScoringConfig class / self.scoring)
        
        Returned scores can then be accessed using NPLinker.all_links or NPLinker.links_for_obj.

        Args:
            objects: An instance of one of the above types, or a list containing any mix of them.
            scoring_method: selects the scoring method to be used (e.g. metcalf, likescore). The value
                of this parameter should be one of the ScoringConfig.<method> members, which also contain the
                configurable parameters for each method.
            datalinks: The DataLinks instance used for scoring. Defaults to None, in which case the
                standard internal DataLinks instance is used. 

        Returns:
            If either a single instance or a list with uniform type is used as input, the return value
            is a list containing the objects for which links were found to satisfy the scoring criteria.
            In the case where a mixed list of objects is used, the return value will instead be a 3 element
            list, containing the individual results for the classes GCF, Spectrum, and MolecularFamily 
            (in that order). 
        """
        # single object, pass it straight through
        if not isinstance(objects, list):
            logger.debug('get_links: single object')
            return self._get_links([objects], datalinks, scoring_method)

        # check if the list contains one or multiple types of object
        type_counts = {NPLinker.OBJ_CLASSES[t]: 0 for t in range(len(NPLinker.OBJ_CLASSES))}
        for o in objects:
            try:
                type_counts[type(o)] += 1
            except KeyError:
                raise Exception('Unsupported type passed to get_links ({})'.format(type(o)))

        if max(type_counts.values()) == len(objects):
            # if the list seems to have a uniform type, again pass straight through
            logger.debug('get_links: uniformly typed list ({})'.format(type(objects[0])))
            return self._get_links(objects, datalinks, scoring_method)

        # finally, if it's a mix of objects, extract them into separate lists and call 
        # the internal method multiple times instead
        obj_lists = {t: [] for t in NPLinker.OBJ_CLASSES}
        for o in objects:
            obj_lists[type(o)].append(o)

        # keep consistent ordering here: GCF, Spectrum, MolecularFamily
        obj_lists = [obj_lists[GCF], obj_lists[Spectrum], obj_lists[MolecularFamily]]

        results = [[] for t in NPLinker.OBJ_CLASSES]
        logger.debug('get_links: multiple lists')
        for i in range(len(NPLinker.OBJ_CLASSES)):
            if len(obj_lists[i]) > 0:
                results[i].append(self._get_links(obj_lists[i], datalinks, scoring_method))

        return results

    def _store_object_links(self, src, dst, score, scoring_method):
        # structure of self._links is:
        # src objects as keys, values are scoring methods
        #   scoring methods as keys, values are dst objects
        #       dst objects as keys, values are link scores
        if src not in self._links:
            self._links[src] = {m: {} for m in SCORING_METHODS}

        if scoring_method.name not in self._links[src]: 
            raise Exception('Should never happen!')

        self._links[src][scoring_method.name][dst] = score

    def _get_links(self, objects, datalinks, scoring_method):
        if len(objects) == 0:
            return []

        if self._linkfinder is None:
            raise Exception('Need to call process_dataset first')

        if scoring_method.name not in SCORING_METHODS:
            raise Exception('unknown scoring method "{}"'.format(scoring_method.name))

        datalinks = datalinks or self._datalinks
        input_type = type(objects[0])

        if scoring_method.name == 'metcalf':
            logger.debug('_get_links for metcalf scoring, standardised={}'.format(scoring_method.standardised))
            if not scoring_method.standardised:
                # get the basic Metcalf scores and carry on
                results = self._linkfinder.get_links(datalinks, objects, scoring_method.name, scoring_method.cutoff)
            else:
                # get the basic Metcalf scores BUT ignore the cutoff value here as it should only be applied to 
                # the final scores, not the original ones
                results = self._linkfinder.get_links(datalinks, objects, scoring_method.name, None)

                # TODO molfam support still missing here (as in data_linking.py)

                # results is a (3, x) array for spec/molfam input, where (1, :) gives src obj ID, (2, :) gives
                # dst obj ID, and (3, :) gives scores
                # for GCFs you get [spec, molfam] with above substructure


                # NOTE: need to use sets here because if the "else" cases are executed can get lots of duplicate
                # object IDs which messes everything up

                # make type1_objects always == spectra. if input is spectra, just use "objects",
                # otherwise build a list using the object IDs in the results array
                type1_objects = objects if input_type == Spectrum else {self.spectra[int(index)] for index in results[0][self.R_DST_ID]}

                # other objs should be "objects" if input is GCF, otherwise create 
                # from the results array as above
                other_objs = objects if input_type == GCF else {self._gcfs[int(index)] for index in results[0][self.R_DST_ID]}

                # apply score update to each spec:gcf pairing and update the corresponding entry in results
                # (the metcalf_expected and metcalf_variance matrices should have been calculated already)
                new_src, new_dst, new_res = [], [], []

                # iterating over spectra
                for i, type1_obj in enumerate(type1_objects):
                    met_strains = len(type1_obj.dataset_strains)
                    # iterating over GCFs
                    for j, other_obj in enumerate(other_objs):
                        gen_strains = len(other_obj.dataset_strains) 

                        # lookup expected + variance values based on strain counts 
                        expected = self._linkfinder.metcalf_expected[met_strains][gen_strains]
                        variance_sqrt = self._linkfinder.metcalf_variance_sqrt[met_strains][gen_strains]

                        k = i if input_type == GCF else j
                        
                        # calculate the final score based on the basic Metcalf score for these two
                        # particular objects
                        final_score = (results[0][self.R_SCORE][k] - expected) / variance_sqrt
                        
                        # finally apply the scoring cutoff and store the result
                        if scoring_method.cutoff is None or (final_score >= scoring_method.cutoff):
                            new_src.append(results[0][self.R_SRC_ID][k])
                            new_dst.append(results[0][self.R_DST_ID][k])
                            new_res.append(final_score)
                    
                # overwrite original "results" with equivalent new data structure
                results = [np.array([new_src, new_dst, new_res])]
                if input_type == GCF:
                    # TODO molfam...
                    results.append(np.zeros((3, 0)))

        elif scoring_method.name == 'likescore':
            results = self._linkfinder.get_links(datalinks, objects, scoring_method.name, scoring_method.cutoff)
        elif scoring_method.name == 'hg':
            results = self._linkfinder.get_links(datalinks, objects, scoring_method.name, scoring_method.prob)
        else:
            raise Exception('Not handled yet')

        scores_found = set()

        if input_type == GCF:
            logger.debug('get_links: result type is GCF, inputs={}, results={}'.format(len(objects), len(results)))
            # for GCF input, results contains two arrays of shape (3, x), 
            # which contain spec-gcf and fam-gcf links respectively 
            result_gcf_spec, result_gcf_fam = results[0], results[1]

            for res, type_ in [(result_gcf_spec, Spectrum), (result_gcf_fam, MolecularFamily)]:
                if res.shape[1] == 0:
                    logger.debug('Found no links for {} GCFs of type {}'.format(len(objects), type_))
                    continue # no results


                # TODO: self._families is never populated atm, so this breaks with MolecularFamily results!
                if type_ == MolecularFamily:
                    logger.warning('FIXME: MolecularFamily support')
                    continue
                logger.debug('Found {} links for {} GCFs of type {}'.format(res.shape[1], len(objects), type_))

                lookup = {gcf.id: objects.index(gcf) for gcf in objects}

                # for each entry in the results (each Spectrum or MolecularFamily)
                for j in range(res.shape[1]):
                    # extract the ID of the object and get the object itself
                    obj_id = int(res[NPLinker.R_DST_ID, j])
                    obj = self._spectra[obj_id] if type_ == Spectrum else self._families[obj_id]

                    # retrieve the GCF object too
                    index = lookup[res[NPLinker.R_SRC_ID, j]]
                    gcf = objects[index]

                    # finally, create the entry in self._links and record that this GCF
                    # has at least one link associated with it
                    self._store_object_links(gcf, obj, res[NPLinker.R_SCORE, j], scoring_method)
                    scores_found.add(gcf)

        else:
            logger.debug('get_links: result type is Spectrum/MolFam, inputs={}, results={}'.format(len(objects), len(results)))
            # for non-GCF input, result is a list containing a single array, shape (3, x)
            # where x is the total number of links found
            results = results[0]
            if results.shape[1] == 0:
                logger.debug('Found no links for {}'.format(objects))
                return []

            logger.debug('Found {} links for {} Spectra/MolFams'.format(results.shape[1], len(objects)))

            lookup = {obj.id: objects.index(obj) for obj in objects}

            # for each entry in the results (each GCF)
            for j in range(results.shape[1]):
                # extract the ID of the GCF and use that to get the object
                gcf_id = int(results[NPLinker.R_DST_ID, j])
                gcf = self._gcfs[gcf_id]

                # retrieve the source object similarly
                index = lookup[results[NPLinker.R_SRC_ID, j]]
                obj = objects[index]

                # finally, create the entry in self._links and record that this Spectrum or
                # MolecularFamily has at least one link associated with it
                self._store_object_links(obj, gcf, results[NPLinker.R_SCORE, j], scoring_method)
                scores_found.add(obj)

        return list(scores_found)

    def get_common_strains(self, objects_a, objects_b, filter_no_shared=True):
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
    def datalinks(self):
        """
        Returns the main DataLinks object used for scoring etc
        """
        return self._datalinks

    @property
    def scoring(self):
        """
        Returns the ScoringConfig object used to manipulate the scoring process
        """
        return self._config.scoring

    @property
    def repro_data(self):
        """
        Returns the dict of reproducibility data
        """
        return self._repro_data

    @property
    def all_links(self):
        """
        Returns the links for all objects and scoring methods
        """
        return self._links

    def links_for_obj(self, obj, scoring_method, sort=True, type_=None):
        # the type_ part is only really useful for GCFs, to pick between MolFam/Spectrum results
        if obj not in self._links:
            logger.warning('No links found for supplied object. Either get_links was not called first or no results were found')
            return None
        links = [(obj, score) for obj, score in self._links[obj][scoring_method.name].items() if (type_ is None or isinstance(obj, type_))]
        if sort:
            links.sort(key=lambda x: x[1], reverse=True)
        return links

if __name__ == "__main__":
    # can set default logging configuration this way...
    LogConfig.setLogLevel(logging.DEBUG)

    # initialise NPLinker from the command-line arguments
    npl = NPLinker(Args().get_args())

    # load the dataset
    if not npl.load_data():
        print('Failed to load the dataset!')
        sys.exit(-1)

    # make sure metcalf scoring is enabled 
    npl.scoring.metcalf.enabled = True

    # generate all datalinks objects etc
    if not npl.process_dataset():
        print('Failed to process dataset')
        sys.exit(-1)

    # pick a GCF to get links for
    test_gcf = npl.gcfs[8]

    # call get_links and set the scoring method to metcalf by referencing the 
    # scoring.metcalf object
    result = npl.get_links(test_gcf, scoring_method=npl.scoring.metcalf)

    # check if any links were found
    if test_gcf not in result:
        print('No links found!')
        sys.exit(0)

    # retrieve the specific links for this object (get_links() supports passing
    # multiple objects as input and so doesn't return them directly). since multiple
    # scoring methods may be in use, need to supply the desired scoring object as a
    # parameter. By default this will return objects of any available type (so Spectrum
    # or MolecularFamily for GCF input), but the type_ parameter can be used to select
    # a single class if needed.
    test_gcf_links = npl.links_for_obj(test_gcf, npl.scoring.metcalf, type_=Spectrum)

    # print the objects (spectra in this example) and their scores
    for obj, score in test_gcf_links:
        print('{} : {}'.format(obj, score))
        # also retrieve common strains. the return value here is a dict indexed
        # by object tuples (either (Spectrum, GCF) or (MolecularFamily, GCF) depending
        # on the inputs), with the values being lists of shared strains (strings)
        common_strains = npl.get_common_strains(test_gcf, obj)
        strain_names = list(common_strains.values())[0]
        if len(strain_names) == 0:
            print('   No shared strains!')
        else:
            print('   {} shared strains: {}'.format(len(strain_names), strain_names))

    print('{} total links found'.format(len(test_gcf_links)))
