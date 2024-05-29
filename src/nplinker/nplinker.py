from __future__ import annotations
import logging
import sys
from os import PathLike
from pprint import pformat
from typing import TYPE_CHECKING
from . import setup_logging
from .arranger import DatasetArranger
from .config import load_config
from .genomics import BGC
from .genomics import GCF
from .loader import NPLINKER_APP_DATA_DIR
from .loader import DatasetLoader
from .metabolomics import MolecularFamily
from .metabolomics import Spectrum
from .pickler import save_pickled_data
from .scoring.link_collection import LinkCollection
from .scoring.metcalf_scoring import MetcalfScoring
from .scoring.methods import ScoringMethod
from .scoring.np_class_scoring import NPClassScoring
from .scoring.rosetta_scoring import RosettaScoring


if TYPE_CHECKING:
    from collections.abc import Sequence
    from .strain import Strain

logger = logging.getLogger(__name__)


class NPLinker:
    """Main class for the NPLinker application."""

    # allowable types for objects to be passed to scoring methods
    OBJ_CLASSES = [Spectrum, MolecularFamily, GCF, BGC]
    # default set of enabled scoring methods
    # TODO: ideally these shouldn't be hardcoded like this
    SCORING_METHODS = {
        MetcalfScoring.NAME: MetcalfScoring,
        RosettaScoring.NAME: RosettaScoring,
        NPClassScoring.NAME: NPClassScoring,
    }

    def __init__(self, config_file: str | PathLike):
        """Initialise an NPLinker instance."""
        self.config = load_config(config_file)

        setup_logging(
            level=self.config.log.level,
            file=self.config.log.get("file", ""),
            use_console=self.config.log.use_console,
        )
        logger.info(
            "Configuration:\n %s", pformat(self.config.as_dict(), width=20, sort_dicts=False)
        )

        self._loader = DatasetLoader()

        self._spectra = []
        self._bgcs = []
        self._gcfs = []
        self._strains = None
        self._metadata = {}
        self._molfams = []
        self._mibig_bgcs = []
        self._chem_classes = None
        self._class_matches = None

        self._bgc_lookup = {}
        self._gcf_lookup = {}
        self._spec_lookup = {}
        self._mf_lookup = {}

        self._scoring_methods = {}
        config_methods = self.config.get("scoring_methods", [])
        for name, method in NPLinker.SCORING_METHODS.items():
            if len(config_methods) == 0 or name in config_methods:
                self._scoring_methods[name] = method
                logger.info(f"Enabled scoring method: {name}")

        self._scoring_methods_setup_complete = {
            name: False for name in self._scoring_methods.keys()
        }

        self._datalinks = None

        self._repro_data = {}
        repro_file = self.config.get("repro_file")
        if repro_file:
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
        self._repro_data["args"] = {}
        for i, arg in enumerate(sys.argv):
            self._repro_data["args"][i] = arg

        # TODO anything else to include here?

        return self._repro_data

    def save_repro_data(self, filename):
        self._collect_repro_data()
        with open(filename, "wb") as repro_file:
            # TODO is pickle the best format to use?
            save_pickled_data(self._repro_data, repro_file)
            logger.info(f"Saving reproducibility data to {filename}")

    @property
    def root_dir(self) -> str:
        """Returns path to the current dataset root directory.

        Returns:
            The path to the dataset root directory currently in use
        """
        return self.config.root_dir

    @property
    def data_dir(self):
        """Returns path to nplinker/data directory (files packaged with the app itself)."""
        return NPLINKER_APP_DATA_DIR

    @property
    def bigscape_cutoff(self):
        """Returns the current BiGSCAPE clustering cutoff value."""
        return self.config.bigscape.cutoff

    def load_data(self):
        """Loads the basic components of a dataset."""
        arranger = DatasetArranger()
        arranger.arrange()
        self._loader.load()

        self._spectra = self._loader.spectra
        self._molfams = self._loader.molfams
        self._bgcs = self._loader.bgcs
        self._gcfs = self._loader.gcfs
        self._mibig_bgcs = self._loader.mibig_bgcs
        self._strains = self._loader.strains
        self._product_types = self._loader.product_types
        self._chem_classes = self._loader.chem_classes
        self._class_matches = self._loader.class_matches

    # TODO CG: refactor this method and update its unit tests
    def get_links(
        self, input_objects: list, scoring_methods: list, and_mode: bool = True
    ) -> LinkCollection:
        """Find links for a set of input objects (BGCs/GCFs/Spectra/MolFams).

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
            and_mode: determines how results from multiple methods are combined.
                This is ignored if a single method is supplied. If multiple methods
                are used and ``and_mode`` is True, the results will only contain
                links found by ALL methods. If False, results will contain links
                found by ANY method.

        Returns:
            An instance of ``nplinker.scoring.methods.LinkCollection``
        """
        if isinstance(input_objects, list) and len(input_objects) == 0:
            raise Exception("input_objects length must be > 0")

        if isinstance(scoring_methods, list) and len(scoring_methods) == 0:
            raise Exception("scoring_methods length must be > 0")

        # for convenience convert a single scoring object into a single entry
        # list
        if not isinstance(scoring_methods, list):
            scoring_methods = [scoring_methods]

        # check if input_objects is a list of lists. if so there should be one
        # entry for each supplied method for it to be a valid parameter
        if isinstance(input_objects[0], list):
            if len(input_objects) != len(scoring_methods):
                raise Exception(
                    "Number of input_objects lists must match number of scoring_methods (found: {}, expected: {})".format(
                        len(input_objects), len(scoring_methods)
                    )
                )

        # TODO check scoring_methods only contains ScoringMethod-derived
        # instances

        # want everything to be in lists of lists
        if not isinstance(input_objects, list) or (
            isinstance(input_objects, list) and not isinstance(input_objects[0], list)
        ):
            input_objects = [input_objects]

        logger.debug(
            "get_links: {} object sets, {} methods".format(len(input_objects), len(scoring_methods))
        )

        # copy the object set if required to make up the numbers
        if len(input_objects) != len(scoring_methods):
            if len(scoring_methods) < len(input_objects):
                raise Exception("Number of scoring methods must be >= number of input object sets")
            elif (len(scoring_methods) > len(input_objects)) and len(input_objects) != 1:
                raise Exception(
                    "Mismatch between number of scoring methods and input objects ({} vs {})".format(
                        len(scoring_methods), len(input_objects)
                    )
                )
            elif len(scoring_methods) > len(input_objects):
                # this is a special case for convenience: pass in 1 set of objects and multiple methods,
                # result is that set is used for all methods
                logger.debug("Duplicating input object set")
                while len(input_objects) < len(scoring_methods):
                    input_objects.append(input_objects[0])
                    logger.debug("Duplicating input object set")

        link_collection = LinkCollection(and_mode)

        for i, method in enumerate(scoring_methods):
            # do any one-off initialisation required by this method
            if not self._scoring_methods_setup_complete[method.name]:
                logger.debug(f"Doing one-time setup for {method.name}")
                self._scoring_methods[method.name].setup(self)
                self._scoring_methods_setup_complete[method.name] = True

            # should construct a dict of {object_with_link: <link_data>}
            # entries
            objects_for_method = input_objects[i]
            logger.debug(
                "Calling scoring method {} on {} objects".format(
                    method.name, len(objects_for_method)
                )
            )
            link_collection = method.get_links(*objects_for_method, link_collection=link_collection)

        if not self._datalinks:
            logger.debug("Creating internal datalinks object")
            self._datalinks = self.scoring_method(MetcalfScoring.NAME).datalinks
            logger.debug("Created internal datalinks object")

        if len(link_collection) == 0:
            logger.debug("No links found or remaining after merging all method results!")

        # populate shared strain info
        logger.debug("Calculating shared strain information...")
        # TODO more efficient version?
        for source, link_data in link_collection.links.items():
            if isinstance(source, BGC):
                logger.debug("Cannot determine shared strains for BGC input!")
                break

            targets = list(filter(lambda x: not isinstance(x, BGC), link_data.keys()))
            if len(targets) > 0:
                if isinstance(source, GCF):
                    shared_strains = self._datalinks.get_common_strains(targets, [source], True)
                    for target, link in link_data.items():
                        if (target, source) in shared_strains:
                            link.shared_strains = shared_strains[(target, source)]
                else:
                    shared_strains = self._datalinks.get_common_strains([source], targets, True)
                    for target, link in link_data.items():
                        if (source, target) in shared_strains:
                            link.shared_strains = shared_strains[(source, target)]

        logger.info("Finished calculating shared strain information")

        logger.info("Final size of link collection is {}".format(len(link_collection)))
        return link_collection

    def get_common_strains(
        self,
        met: Sequence[Spectrum] | Sequence[MolecularFamily],
        gcfs: Sequence[GCF],
        filter_no_shared: bool = True,
    ) -> dict[tuple[Spectrum | MolecularFamily, GCF], list[Strain]]:
        """Get common strains between given spectra/molecular families and GCFs.

        Args:
            met:
                A list of Spectrum or MolecularFamily objects.
            gcfs: A list of GCF objects.
            filter_no_shared: If True, the pairs of spectrum/mf and GCF
                without common strains will be removed from the returned dict;

        Returns:
            A dict where the keys are tuples of (Spectrum/MolecularFamily, GCF)
            and values are a list of shared Strain objects.
        """
        if not self._datalinks:
            self._datalinks = self.scoring_method(MetcalfScoring.NAME).datalinks
        common_strains = self._datalinks.get_common_strains(met, gcfs, filter_no_shared)
        return common_strains

    def has_bgc(self, bgc_id):
        """Returns True if BGC ``bgc_id`` exists in the dataset."""
        return bgc_id in self._bgc_lookup

    def lookup_bgc(self, bgc_id):
        """If BGC ``bgc_id`` exists, return it. Otherwise return None."""
        return self._bgc_lookup.get(bgc_id, None)

    def lookup_gcf(self, gcf_id):
        """If GCF ``gcf_id`` exists, return it. Otherwise return None."""
        return self._gcf_lookup.get(gcf_id, None)

    def lookup_spectrum(self, spectrum_id):
        """If Spectrum ``name`` exists, return it. Otherwise return None."""
        return self._spec_lookup.get(spectrum_id, None)

    def lookup_mf(self, mf_id):
        """If MolecularFamily `family_id` exists, return it. Otherwise return None."""
        return self._mf_lookup.get(mf_id, None)

    @property
    def strains(self):
        """Returns a list of all the strains in the dataset."""
        return self._strains

    @property
    def bgcs(self):
        """Returns a list of all the BGCs in the dataset."""
        return self._bgcs

    @property
    def gcfs(self):
        """Returns a list of all the GCFs in the dataset."""
        return self._gcfs

    @property
    def spectra(self):
        """Returns a list of all the Spectra in the dataset."""
        return self._spectra

    @property
    def molfams(self):
        """Returns a list of all the MolecularFamilies in the dataset."""
        return self._molfams

    @property
    def metadata(self):
        return self._metadata

    @property
    def mibig_bgcs(self):
        """Get a list of all the MIBiG BGCs in the dataset."""
        return self._mibig_bgcs

    @property
    def product_types(self):
        """Returns a list of the available BiGSCAPE product types in current dataset."""
        return self._product_types

    @property
    def repro_data(self):
        """Returns the dict containing reproducibility data."""
        return self._repro_data

    @property
    def scoring_methods(self):
        """Returns a list of available scoring method names."""
        return list(self._scoring_methods.keys())

    @property
    def chem_classes(self):
        """Returns loaded ChemClassPredictions with the class predictions."""
        return self._chem_classes

    @property
    def class_matches(self):
        """ClassMatches with the matched classes and scoring tables from MIBiG."""
        return self._class_matches

    def scoring_method(self, name: str) -> ScoringMethod | None:
        """Return an instance of a scoring method.

        Args:
            name: the name of the method (see :func:`scoring_methods`)

        Returns:
            An instance of the named scoring method class, or None if the name is invalid
        """
        if name not in self._scoring_methods_setup_complete:
            return None

        if not self._scoring_methods_setup_complete[name]:
            self._scoring_methods[name].setup(self)
            self._scoring_methods_setup_complete[name] = True

        return self._scoring_methods.get(name, None)(self)
