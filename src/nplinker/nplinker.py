from __future__ import annotations
import copy
import logging
import sys
from typing import TYPE_CHECKING
from .config import Config
from .genomics import BGC
from .genomics import GCF
from .loader import NPLINKER_APP_DATA_DIR
from .loader import DatasetLoader
from .logconfig import LogConfig
from .metabolomics.molecular_family import MolecularFamily
from .metabolomics.spectrum import Spectrum
from .pickler import save_pickled_data
from .scoring.link_collection import LinkCollection
from .scoring.metcalf_scoring import MetcalfScoring
from .scoring.np_class_scoring import NPClassScoring
from .scoring.rosetta_scoring import RosettaScoring
from .strain_collection import StrainCollection


if TYPE_CHECKING:
    from collections.abc import Sequence
    from .strains import Strain

logger = LogConfig.getLogger(__name__)


class NPLinker():

    # allowable types for objects to be passed to scoring methods
    OBJ_CLASSES = [Spectrum, MolecularFamily, GCF, BGC]
    # default set of enabled scoring methods
    # TODO: ideally these shouldn't be hardcoded like this
    SCORING_METHODS = {
        MetcalfScoring.NAME: MetcalfScoring,
        RosettaScoring.NAME: RosettaScoring,
        NPClassScoring.NAME: NPClassScoring
    }

    def __init__(self, userconfig=None):
        """Initialise an NPLinker instance.

        NPLinker instances can be configured in multiple ways, in ascending order of priority:

        1. A global user-level default configuration file in TOML format, found in the directory XDG_CONFIG_HOME/nplinker/nplinker.toml
        2. A local TOML configuration file
        3. Command-line arguments / supplying a manually constructed parameter dictionary

        The global user-level configuration file will be created automatically the first time
        an NPLinker instance is initialised if it doesn't already exist. The default file contains
        sensible default values for each setting and is intended to be copied and edited to
        produce dataset-specific configuration files, which will then override any parameters shared
        with the user-level file. To load such a file, simply set the "userconfig" parameter to a
        string containing the filename.

        It's also possible to selectively override configuration file parameters by
        supplying command-line arguments (if running nplinker.py as a script), or by passing
        a dict with a structure corresponding to the configuration file format to this method.

        Some examples may make the various possible combinations a bit clearer::

            # simplest option: load a local configuration file
            > npl = NPLinker('myconfig.toml')

            # the same thing but running as a script
            > python -m nplinker.nplinker --config "myconfig.toml"

            # use the defaults from the user-level config while modifying the root path
            # to load the dataset from (this is the minimum you would need to change in the
            # default config file)
            > npl = NPLinker({'dataset': {'root': '/path/to/dataset'}})

            # the same thing running NPLinker as a script
            > python nplinker.py --dataset.root /path/to/dataset

        Args:
            userconfig: supplies user-defined configuration data. May take one of 3 types:

                - str: treat as filename of a local configuration file to load
                        (overriding the defaults)
                - dict: contents will be used to override values in the dict generated
                        from loading the configuration file(s)
        """

        # if userconfig is a string => create a dict with 'config' key and string as filename
        # if userconfig is a dict => pass it to Config() directly
        if isinstance(userconfig, str):
            userconfig = {'config': userconfig}
        elif not isinstance(userconfig, dict):
            raise Exception(
                'Invalid type for userconfig (should be None/str/dict, found "{}")'
                .format(type(userconfig)))

        self._config = Config(userconfig)

        # configure logging based on the supplied config params
        LogConfig.setLogLevelStr(self._config.config['loglevel'])
        logfile = self._config.config['logfile']
        if len(logfile) > 0:
            logfile_dest = logging.FileHandler(logfile)
            # if we want to log to stdout plus logfile, add the new destination
            if self._config.config.get('log_to_stdout',
                                       True):  # default to True
                LogConfig.addLogDestination(logfile_dest)
            else:
                # otherwise overwrite the default stdout destination
                LogConfig.setLogDestination(logfile_dest)

        # the DatasetLoader takes care of figuring out the locations of all the relevant files/folders
        # and will show error/warning if any missing (depends if optional or
        # not)
        self._loader = DatasetLoader(self._config.config)

        self._spectra = []
        self._bgcs = []
        self._gcfs = []
        self._strains = None
        self._metadata = {}
        self._molfams = []
        self._mibig_bgc_dict = {}
        self._chem_classes = None
        self._class_matches = None

        self._bgc_lookup = {}
        self._gcf_lookup = {}
        self._spec_lookup = {}
        self._mf_lookup = {}

        self._scoring_methods = {}
        config_methods = self._config.config.get('scoring_methods', [])
        for name, method in NPLinker.SCORING_METHODS.items():
            if len(config_methods) == 0 or name in config_methods:
                self._scoring_methods[name] = method
                logger.debug(f'Enabled scoring method: {name}')

        self._scoring_methods_setup_complete = {
            name: False
            for name in self._scoring_methods.keys()
        }

        self._datalinks = None

        self._repro_data = {}
        repro_file = self._config.config['repro_file']
        if len(repro_file) > 0:
            self.save_repro_data(repro_file)

    def _collect_repro_data(self):
        """Creates a dict containing data to aid reproducible runs of nplinker.

        This method creates a dict containing various bits of information about
        the current execution of nplinker. This data will typically be saved to
        a file in order to aid reproducibility using :func:`save_repro_data`.

        TODO describe contents

        Returns:
            A dict containing the information described above
        """

        self._repro_data = {}
        # TODO best way to embed latest git commit hash? probably with a packaging script...
        # TODO versions of all Python dependencies used (just rely on
        # Pipfile.lock here?)

        # insert command line arguments
        self._repro_data['args'] = {}
        for i, arg in enumerate(sys.argv):
            self._repro_data['args'][i] = arg

        # TODO anything else to include here?

        return self._repro_data

    def save_repro_data(self, filename):
        self._collect_repro_data()
        with open(filename, 'wb') as repro_file:
            # TODO is pickle the best format to use?
            save_pickled_data(self._repro_data, repro_file)
            logger.info(f'Saving reproducibility data to {filename}')

    @property
    def config(self):
        """Returns a copy of the data parsed from the configuration file as a dict

        Returns:
                dict: configuration file parameters as a nested dict
        """
        return copy.deepcopy(self._config.config)

    @property
    def root_dir(self):
        """Returns path to the current dataset root directory

        Returns:
                str: the path to the dataset root directory currently in use
        """
        return self._loader._root

    @property
    def dataset_id(self):
        """Returns dataset "ID".

        For local datasets this will just be the last component of the directory path,
        e.g. /path/to/my_dataset would produce an ID of "my_dataset".

        For datasets loaded from the Paired Omics platform the ID will be the platform
        project ID, e.g. "MSV000079284"

        Returns:
            str: the dataset ID
        """
        return self._loader.dataset_id

    @property
    def data_dir(self):
        """Returns path to nplinker/data directory (files packaged with the app itself)"""
        return NPLINKER_APP_DATA_DIR

    @property
    def gnps_params(self):
        """Returns a dict containing data from GNPS params.xml (if available).

        Returns:
            dict: GNPS parameters, or an empty dict if none exist in the dataset
        """
        return self._loader.gnps_params

    @property
    def dataset_description(self):
        """Returns dataset description.

        If nplinker finds a 'description.txt' file in the root directory of the
        dataset, the content will be parsed and made available through this property.

        Returns:
            str: the content of description.txt or '<no description>'
        """
        return self._loader.description_text

    @property
    def bigscape_cutoff(self):
        """Returns the current BiGSCAPE clustering cutoff value"""
        return self._loader._bigscape_cutoff

    def load_data(self, new_bigscape_cutoff=None):
        """Loads the basic components of a dataset.

        This method is responsible for loading the various pieces of the supplied dataset into
        memory and doing any initial parsing/object manipulation required. After it completes,
        applications can access the lists of GCFs, Spectra, MolecularFamilies and strains
        using the corresponding properties of the NPLinker class.

        Returns:
            bool: True if successful, False otherwise
        """
        logger.debug('load_data(new_bigscape_cutoff=%s)', new_bigscape_cutoff)
        if new_bigscape_cutoff is None:
            self._loader.validate()
            if not self._loader.load():
                return False
        else:
            # CG: only reload genomics data when changing bigscape cutoff
            # TODO: this part should be removed, reload everything if bigscape data changes.
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
        self._chem_classes = self._loader.chem_classes
        self._class_matches = self._loader.class_matches

        logger.debug('Generating lookup tables: genomics')
        self._bgc_lookup = {bgc.bgc_id: bgc for bgc in self._bgcs}
        self._gcf_lookup = {gcf.gcf_id: gcf for gcf in self._gcfs}

        # don't need to do these two if cutoff changed (indicating genomics data
        # was reloaded but not metabolomics)
        if new_bigscape_cutoff is None:
            logger.debug('Generating lookup tables: metabolomics')
            self._spec_lookup = {
                spec.spectrum_id: spec
                for spec in self._spectra
            }
            self._mf_lookup = {mf.family_id: mf for mf in self._molfams}

        logger.debug('load_data: completed')
        return True

    # TODO CG: refactor this method and update its unit tests
    def get_links(self, input_objects, scoring_methods, and_mode=True):
        """Find links for a set of input objects (BGCs/GCFs/Spectra/MolFams)

        The input objects can be any mix of the following NPLinker types:

            - BGC
            - GCF
            - Spectrum
            - MolecularFamily

        TODO longer description here

        Args:
            input_objects: objects to be passed to the scoring method(s).
                This may be either a flat list of a uniform type (one of the 4
                types above), or a list of such lists
            scoring_methods: a list of one or more scoring methods to use
            and_mode (bool): determines how results from multiple methods are combined.
                This is ignored if a single method is supplied. If multiple methods
                are used and ``and_mode`` is True, the results will only contain
                links found by ALL methods. If False, results will contain links
                found by ANY method.

        Returns:
            An instance of ``nplinker.scoring.methods.LinkCollection``
        """
        if isinstance(input_objects, list) and len(input_objects) == 0:
            raise Exception('input_objects length must be > 0')

        if isinstance(scoring_methods, list) and len(scoring_methods) == 0:
            raise Exception('scoring_methods length must be > 0')

        # for convenience convert a single scoring object into a single entry
        # list
        if not isinstance(scoring_methods, list):
            scoring_methods = [scoring_methods]

        # check if input_objects is a list of lists. if so there should be one
        # entry for each supplied method for it to be a valid parameter
        if isinstance(input_objects[0], list):
            if len(input_objects) != len(scoring_methods):
                raise Exception(
                    'Number of input_objects lists must match number of scoring_methods (found: {}, expected: {})'
                    .format(len(input_objects), len(scoring_methods)))

        # TODO check scoring_methods only contains ScoringMethod-derived
        # instances

        # want everything to be in lists of lists
        if not isinstance(input_objects, list) or (
                isinstance(input_objects, list)
                and not isinstance(input_objects[0], list)):
            input_objects = [input_objects]

        logger.debug('get_links: {} object sets, {} methods'.format(
            len(input_objects), len(scoring_methods)))

        # copy the object set if required to make up the numbers
        if len(input_objects) != len(scoring_methods):
            if len(scoring_methods) < len(input_objects):
                raise Exception(
                    'Number of scoring methods must be >= number of input object sets'
                )
            elif (len(scoring_methods) >
                  len(input_objects)) and len(input_objects) != 1:
                raise Exception(
                    'Mismatch between number of scoring methods and input objects ({} vs {})'
                    .format(len(scoring_methods), len(input_objects)))
            elif len(scoring_methods) > len(input_objects):
                # this is a special case for convenience: pass in 1 set of objects and multiple methods,
                # result is that set is used for all methods
                logger.debug('Duplicating input object set')
                while len(input_objects) < len(scoring_methods):
                    input_objects.append(input_objects[0])
                    logger.debug('Duplicating input object set')

        link_collection = LinkCollection(and_mode)

        for i, method in enumerate(scoring_methods):
            # do any one-off initialisation required by this method
            if not self._scoring_methods_setup_complete[method.name]:
                logger.debug(f'Doing one-time setup for {method.name}')
                self._scoring_methods[method.name].setup(self)
                self._scoring_methods_setup_complete[method.name] = True

            # should construct a dict of {object_with_link: <link_data>}
            # entries
            objects_for_method = input_objects[i]
            logger.debug('Calling scoring method {} on {} objects'.format(
                method.name, len(objects_for_method)))
            link_collection = method.get_links(*objects_for_method,
                                               link_collection=link_collection)

        if not self._datalinks:
            logger.debug('Creating internal datalinks object')
            self._datalinks = self.scoring_method(
                MetcalfScoring.NAME).datalinks
            logger.debug('Created internal datalinks object')

        if len(link_collection) == 0:
            logger.debug(
                'No links found or remaining after merging all method results!'
            )

        # populate shared strain info
        logger.debug('Calculating shared strain information...')
        # TODO more efficient version?
        for source, link_data in link_collection.links.items():
            if isinstance(source, BGC):
                logger.debug('Cannot determine shared strains for BGC input!')
                break

            targets = list(
                filter(lambda x: not isinstance(x, BGC), link_data.keys()))
            if len(targets) > 0:
                if isinstance(source, GCF):
                    shared_strains = self._datalinks.get_common_strains(
                        targets, [source], True)
                    for target, link in link_data.items():
                        if (target, source) in shared_strains:
                            link.shared_strains = shared_strains[(target,
                                                                  source)]
                else:
                    shared_strains = self._datalinks.get_common_strains(
                        [source], targets, True)
                    for target, link in link_data.items():
                        if (source, target) in shared_strains:
                            link.shared_strains = shared_strains[(source,
                                                                  target)]

        logger.debug('Finished calculating shared strain information')

        logger.debug('Final size of link collection is {}'.format(
            len(link_collection)))
        return link_collection

    def get_common_strains(
        self,
        met: Sequence[Spectrum] | Sequence[MolecularFamily],
        gcfs: Sequence[GCF],
        filter_no_shared: bool = True
    ) -> dict[tuple[Spectrum | MolecularFamily, GCF], list[Strain]]:
        """Get common strains between given spectra/molecular families and GCFs.

        Note that SingletonFamily objects are excluded from given molecular families.

        Args:
            met(Sequence[Spectrum] | Sequence[MolecularFamily]):
                A list of Spectrum or MolecularFamily objects.
            gcfs(Sequence[GCF]): A list of GCF objects.
            filter_no_shared(bool): If True, the pairs of spectrum/mf and GCF
                without common strains will be removed from the returned dict;

        Returns:
            dict: A dict where the keys are tuples of (Spectrum/MolecularFamily, GCF)
            and values are a list of shared Strain objects.
        """
        if not self._datalinks:
            self._datalinks = self.scoring_method(
                MetcalfScoring.NAME).datalinks
        common_strains = self._datalinks.get_common_strains(
            met, gcfs, filter_no_shared)
        return common_strains

    def has_bgc(self, bgc_id):
        """Returns True if BGC ``bgc_id`` exists in the dataset"""
        return bgc_id in self._bgc_lookup

    def lookup_bgc(self, bgc_id):
        """If BGC ``bgc_id`` exists, return it. Otherwise return None"""
        return self._bgc_lookup.get(bgc_id, None)

    def lookup_gcf(self, gcf_id):
        """If GCF ``gcf_id`` exists, return it. Otherwise return None"""
        return self._gcf_lookup.get(gcf_id, None)

    def lookup_spectrum(self, spectrum_id):
        """If Spectrum ``name`` exists, return it. Otherwise return None"""
        return self._spec_lookup.get(spectrum_id, None)

    def lookup_mf(self, mf_id):
        """If MolecularFamily `family_id` exists, return it. Otherwise return None"""
        return self._mf_lookup.get(mf_id, None)

    @property
    def strains(self):
        """Returns a list of all the strains in the dataset"""
        return self._strains

    @property
    def bgcs(self):
        """Returns a list of all the BGCs in the dataset"""
        return self._bgcs

    @property
    def gcfs(self):
        """Returns a list of all the GCFs in the dataset"""
        return self._gcfs

    @property
    def spectra(self):
        """Returns a list of all the Spectra in the dataset"""
        return self._spectra

    @property
    def molfams(self):
        """Returns a list of all the MolecularFamilies in the dataset"""
        return self._molfams

    @property
    def metadata(self):
        return self._metadata

    @property
    def mibig_bgc_dict(self):
        return self._mibig_bgc_dict

    @property
    def product_types(self):
        """Returns a list of the available BiGSCAPE product types in current dataset"""
        return self._product_types

    @property
    def repro_data(self):
        """Returns the dict containing reproducibility data"""
        return self._repro_data

    @property
    def scoring_methods(self):
        """Returns a list of available scoring method names"""
        return list(self._scoring_methods.keys())

    @property
    def chem_classes(self):
        """Returns loaded ChemClassPredictions with the class predictions"""
        return self._chem_classes

    @property
    def class_matches(self):
        """ClassMatches with the matched classes and scoring tables from MIBiG
        """
        return self._class_matches

    def scoring_method(self, name):
        """Return an instance of a scoring method.

        Args:
            name (str): the name of the method (see :func:`scoring_methods`)

        Returns:
            An instance of the named scoring method class, or None if the name is invalid
        """
        if name not in self._scoring_methods_setup_complete:
            return None

        if not self._scoring_methods_setup_complete[name]:
            self._scoring_methods[name].setup(self)
            self._scoring_methods_setup_complete[name] = True

        return self._scoring_methods.get(name, None)(self)
