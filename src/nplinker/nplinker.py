from __future__ import annotations
import logging
from os import PathLike
from pprint import pformat
from nplinker.strain.strain_collection import StrainCollection
from . import setup_logging
from .arranger import DatasetArranger
from .config import load_config
from .defaults import OUTPUT_DIRNAME
from .genomics import BGC
from .genomics import GCF
from .loader import DatasetLoader
from .metabolomics import MolecularFamily
from .metabolomics import Spectrum
from .scoring.abc import ScoringBase
from .scoring.metcalf_scoring import MetcalfScoring
from .scoring.np_class_scoring import NPClassScoring
from .scoring.rosetta_scoring import RosettaScoring


logger = logging.getLogger(__name__)


class NPLinker:
    """Main class for the NPLinker application."""

    # allowable types for objects to be passed to scoring methods
    OBJ_CLASSES = [Spectrum, MolecularFamily, GCF, BGC]
    # default set of enabled scoring methods
    # TODO: ideally these shouldn't be hardcoded like this
    SCORING_METHODS = {
        MetcalfScoring.name: MetcalfScoring,
        RosettaScoring.name: RosettaScoring,
        NPClassScoring.name: NPClassScoring,
    }

    def __init__(self, config_file: str | PathLike):
        """Initialise an NPLinker instance.

        Args:
            config_file: Path to the configuration file to use.
        """
        self.config = load_config(config_file)

        setup_logging(
            level=self.config.log.level,
            file=self.config.log.get("file", ""),
            use_console=self.config.log.use_console,
        )
        logger.info(
            "Configuration:\n %s", pformat(self.config.as_dict(), width=20, sort_dicts=False)
        )

        self.output_dir = self.config.root_dir / OUTPUT_DIRNAME
        self.output_dir.mkdir(exist_ok=True)

        # data containers that will be populated by the `load_data` method
        self._bgc_dict = {}
        self._gcf_dict = {}
        self._spec_dict = {}
        self._mf_dict = {}
        self._mibig_bgcs = []
        self._strains = StrainCollection()
        self._product_types = []
        self._chem_classes = None
        self._class_matches = None

        self._scoring_methods = {}
        config_methods = self.config.get("scoring_methods", [])
        for name, method in NPLinker.SCORING_METHODS.items():
            if len(config_methods) == 0 or name in config_methods:
                self._scoring_methods[name] = method
                logger.info(f"Enabled scoring method: {name}")

        self._scoring_methods_setup_complete = {
            name: False for name in self._scoring_methods.keys()
        }

    @property
    def root_dir(self) -> str:
        """Returns path to the current dataset root directory.

        Returns:
            The path to the dataset root directory currently in use
        """
        return self.config.root_dir

    @property
    def bigscape_cutoff(self):
        """Returns the current BiGSCAPE clustering cutoff value."""
        return self.config.bigscape.cutoff

    def load_data(self):
        """Load all data from local files into memory.

        This method is a convenience function that calls the `DatasetArranger` and `DatasetLoader`
        classes to load all data from the local filesystem into memory. The loaded data is then
        stored in various private data containers for easy access.
        """
        arranger = DatasetArranger(self.config)
        arranger.arrange()
        loader = DatasetLoader(self.config)
        loader.load()

        self._bgc_dict = {bgc.id: bgc for bgc in loader.bgcs}
        self._gcf_dict = {gcf.id: gcf for gcf in loader.gcfs}
        self._spec_dict = {spec.id: spec for spec in loader.spectra}
        self._mf_dict = {mf.id: mf for mf in loader.mfs}

        self._mibig_bgcs = loader.mibig_bgcs
        self._strains = loader.strains
        self._product_types = loader.product_types
        self._chem_classes = loader.chem_classes
        self._class_matches = loader.class_matches

    # TODO CG: refactor this method and update its unit tests
    def get_links(
        self, input_objects: list, scoring_methods: list, and_mode: bool = True
    ) -> LinkCollection:
        """Find links for a set of input objects (BGCs/GCFs/Spectra/mfs).

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

        if len(link_collection) == 0:
            logger.debug("No links found or remaining after merging all method results!")

        logger.info("Final size of link collection is {}".format(len(link_collection)))
        return link_collection

    def lookup_bgc(self, id: str) -> BGC | None:
        """Get the BGC object with the given ID.

        Args:
            id: the ID of the BGC to look up.

        Returns:
            The BGC object with the given ID, or None if no such object exists.
        """
        return self._bgc_dict.get(id, None)

    def lookup_gcf(self, id: str) -> GCF | None:
        """Get the GCF object with the given ID.

        Args:
            id: the ID of the GCF to look up.

        Returns:
            The GCF object with the given ID, or None if no such object exists.
        """
        return self._gcf_dict.get(id, None)

    def lookup_spectrum(self, id: str) -> Spectrum | None:
        """Get the Spectrum object with the given ID.

        Args:
            id: the ID of the Spectrum to look up.

        Returns:
            The Spectrum object with the given ID, or None if no such object exists.
        """
        return self._spec_dict.get(id, None)

    def lookup_mf(self, id: str) -> MolecularFamily | None:
        """Get the MolecularFamily object with the given ID.

        Args:
            id: the ID of the MolecularFamily to look up.

        Returns:
            The MolecularFamily object with the given ID, or None if no such object exists.
        """
        return self._mf_dict.get(id, None)

    @property
    def strains(self) -> StrainCollection:
        """Get all Strain objects."""
        return self._strains

    @property
    def bgcs(self) -> list[BGC]:
        """Get all BGC objects."""
        return self._bgc_dict.values()

    @property
    def gcfs(self) -> list[GCF]:
        """Get all GCF objects."""
        return self._gcf_dict.values()

    @property
    def spectra(self) -> list[Spectrum]:
        """Get all Spectrum objects."""
        return self._spec_dict.values()

    @property
    def mfs(self) -> list[MolecularFamily]:
        """Get all MolecularFamily objects."""
        return self._mf_dict.values()

    @property
    def mibig_bgcs(self) -> list[BGC]:
        """Get all MiBIG BGC objects."""
        return self._mibig_bgcs

    @property
    def product_types(self) -> list[str]:
        """Get all BiGSCAPE product types."""
        return self._product_types

    @property
    def chem_classes(self):
        """Returns loaded ChemClassPredictions with the class predictions."""
        return self._chem_classes

    @property
    def class_matches(self):
        """ClassMatches with the matched classes and scoring tables from MIBiG."""
        return self._class_matches

    @property
    def scoring_methods(self):
        """Returns a list of available scoring method names."""
        return list(self._scoring_methods.keys())

    def scoring_method(self, name: str) -> ScoringBase | None:
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
