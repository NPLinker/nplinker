from __future__ import annotations
import logging
import pickle
from collections.abc import Sequence
from os import PathLike
from pprint import pformat
from typing import Any
from typing import overload
from dynaconf import Dynaconf
from .arranger import DatasetArranger
from .config import load_config
from .defaults import OUTPUT_DIRNAME
from .genomics import BGC
from .genomics import GCF
from .loader import DatasetLoader
from .logger import setup_logging
from .metabolomics import MolecularFamily
from .metabolomics import Spectrum
from .scoring.link_graph import LinkGraph
from .scoring.metcalf_scoring import MetcalfScoring
from .strain import StrainCollection


logger = logging.getLogger(__name__)


class NPLinker:
    """The central class of NPLinker application.

    Attributes:
        config: The configuration object for the current NPLinker application.
        root_dir: The path to the root directory of the current NPLinker application.
        output_dir: The path to the output directory of the current NPLinker application.
        bgcs: A list of all BGC objects.
        gcfs: A list of all GCF objects.
        spectra: A list of all Spectrum objects.
        mfs: A list of all MolecularFamily objects.
        mibig_bgcs: A list of all MiBIG BGC objects.
        strains: A StrainCollection object containing all Strain objects.
        product_types: A list of all BiGSCAPE product types.
        scoring_methods: A list of all valid scoring methods.
    """

    # Valid scoring methods
    _valid_scoring_methods = {
        MetcalfScoring.name: MetcalfScoring,
        # RosettaScoring.name: RosettaScoring, # To be refactored
        # NPClassScoring.name: NPClassScoring, # To be refactored
    }

    def __init__(self, config_file: str | PathLike):
        """Initialise an NPLinker instance.

        Args:
            config_file: Path to the configuration file to use.


        Examples:
            Starting the NPLinker application:
            >>> from nplinker import NPLinker
            >>> npl = NPLinker("path/to/config.toml")

            Loading data from files to python objects:
            >>> npl.load_data()

            Checking the number of GCF objects:
            >>> len(npl.gcfs)

            Getting the links for all GCF objects using the Metcalf scoring method, and the result
            is stored in a [LinkGraph][nplinker.scoring.LinkGraph] object:
            >>> lg = npl.get_links(npl.gcfs, "metcalf")

            Getting the link data between two objects:
            >>> link_data = lg.get_link_data(npl.gcfs[0], npl.spectra[0])
            {"metcalf": Score("metcalf", 1.0, {"cutoff": 0, "standardised": False})}

            Saving the data to a pickle file:
            >>> npl.save_data("path/to/output.pkl", lg)
        """
        # Load the configuration file
        self.config: Dynaconf = load_config(config_file)

        # Setup logging for the application
        setup_logging(
            level=self.config.log.level,
            file=self.config.log.get("file", ""),
            use_console=self.config.log.use_console,
        )
        logger.info(
            "Configuration:\n %s", pformat(self.config.as_dict(), width=20, sort_dicts=False)
        )

        # Setup the output directory
        self._output_dir = self.config.root_dir / OUTPUT_DIRNAME
        self._output_dir.mkdir(exist_ok=True)

        # Initialise data containers that will be populated by the `load_data` method
        self._bgc_dict: dict[str, BGC] = {}
        self._gcf_dict: dict[str, GCF] = {}
        self._spec_dict: dict[str, Spectrum] = {}
        self._mf_dict: dict[str, MolecularFamily] = {}
        self._mibig_bgcs: list[BGC] = []
        self._strains: StrainCollection = StrainCollection()
        self._product_types: list = []
        self._chem_classes = None  # TODO: to be refactored
        self._class_matches = None  # TODO: to be refactored

        # Flags to keep track of whether the scoring methods have been set up
        self._scoring_methods_setup_done = {name: False for name in self._valid_scoring_methods}

    @property
    def root_dir(self) -> str:
        """Get the path to the root directory of the current NPLinker instance."""
        return str(self.config.root_dir)

    @property
    def output_dir(self) -> str:
        """Get the path to the output directory of the current NPLinker instance."""
        return str(self._output_dir)

    @property
    def bgcs(self) -> list[BGC]:
        """Get all BGC objects."""
        return list(self._bgc_dict.values())

    @property
    def gcfs(self) -> list[GCF]:
        """Get all GCF objects."""
        return list(self._gcf_dict.values())

    @property
    def spectra(self) -> list[Spectrum]:
        """Get all Spectrum objects."""
        return list(self._spec_dict.values())

    @property
    def mfs(self) -> list[MolecularFamily]:
        """Get all MolecularFamily objects."""
        return list(self._mf_dict.values())

    @property
    def mibig_bgcs(self) -> list[BGC]:
        """Get all MiBIG BGC objects."""
        return self._mibig_bgcs

    @property
    def strains(self) -> StrainCollection:
        """Get all Strain objects."""
        return self._strains

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
    def scoring_methods(self) -> list[str]:
        """Get names of all valid scoring methods."""
        return list(self._valid_scoring_methods.keys())

    def load_data(self):
        """Load all data from files into memory.

        This method is a convenience function that calls the
        [`DatasetArranger`][nplinker.arranger.DatasetArranger] class to arrange data files
        (download, generate and/or validate data) in the [correct directory structure][working-directory-structure],
        and then calls the [`DatasetLoader`][nplinker.loader.DatasetLoader] class to load all data
        from the files into memory.

        The loaded data is stored in various data containers for easy access, e.g.
        [`self.bgcs`][nplinker.NPLinker.bgcs] for all BGC objects,
        [`self.strains`][nplinker.NPLinker.strains] for all Strain objects, etc.
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

    @overload
    def get_links(
        self, objects: Sequence[BGC], scoring_method: str, **scoring_params: Any
    ) -> LinkGraph: ...
    @overload
    def get_links(
        self, objects: Sequence[GCF], scoring_method: str, **scoring_params: Any
    ) -> LinkGraph: ...
    @overload
    def get_links(
        self, objects: Sequence[Spectrum], scoring_method: str, **scoring_params: Any
    ) -> LinkGraph: ...
    @overload
    def get_links(
        self, objects: Sequence[MolecularFamily], scoring_method: str, **scoring_params: Any
    ) -> LinkGraph: ...

    def get_links(
        self,
        objects: Sequence[BGC] | Sequence[GCF] | Sequence[Spectrum] | Sequence[MolecularFamily],
        scoring_method: str,
        **scoring_params: Any,
    ) -> LinkGraph:
        """Get links for the given objects using the specified scoring method and parameters.

        Args:
            objects: A sequence of objects to get links for. The objects must be of the same
                type, i.e. `BGC`, `GCF`, `Spectrum` or `MolecularFamily` type.
                !!! Warning
                    For scoring method `metcalf`, the `BGC` objects are not supported.
            scoring_method: The scoring method to use. Must be one of the valid scoring methods
                [`self.scoring_methods`][nplinker.NPLinker.scoring_methods], such as `metcalf`.
            scoring_params: Parameters to pass to the scoring method. If not given, the default
                parameters of the specified scoring method will be used.

                Check the `get_links` method of the scoring method class for the available
                parameters and their default values.

                | Scoring Method | Scoring Parameters |
                | -------------- | ------------------ |
                | `metcalf` | [`cutoff`, `standardised`][nplinker.scoring.MetcalfScoring.get_links] |

        Returns:
            A LinkGraph object containing the links for the given objects.

        Raises:
            ValueError: If input objects are empty or if the scoring method is invalid.
            TypeError: If the input objects are not of the same type or if the object type is invalid.

        Examples:
            Using default scoring parameters:
            >>> lg = npl.get_links(npl.gcfs, "metcalf")

            Scoring parameters provided:
            >>> lg = npl.get_links(npl.gcfs, "metcalf", cutoff=0.5, standardised=True)
        """
        # Validate objects
        if len(objects) == 0:
            raise ValueError("No objects provided to get links for")
        # check if all objects are of the same type
        types = {type(i) for i in objects}
        if len(types) > 1:
            raise TypeError("Input objects must be of the same type.")
        # check if the object type is valid
        obj_type = next(iter(types))
        if obj_type not in (BGC, GCF, Spectrum, MolecularFamily):
            raise TypeError(
                f"Invalid type {obj_type}. Input objects must be BGC, GCF, Spectrum or MolecularFamily objects."
            )

        # Validate scoring method
        if scoring_method not in self._valid_scoring_methods:
            raise ValueError(f"Invalid scoring method {scoring_method}.")

        # Check if the scoring method has been set up
        if not self._scoring_methods_setup_done[scoring_method]:
            self._valid_scoring_methods[scoring_method].setup(self)
            self._scoring_methods_setup_done[scoring_method] = True

        # Initialise the scoring method
        scoring = self._valid_scoring_methods[scoring_method]()

        return scoring.get_links(*objects, **scoring_params)

    def lookup_bgc(self, id: str) -> BGC | None:
        """Get the BGC object with the given ID.

        Args:
            id: the ID of the BGC to look up.

        Returns:
            The BGC object with the given ID, or None if no such object exists.

        Examples:
            >>> bgc = npl.lookup_bgc("BGC000001")
            >>> bgc
            BGC(id="BGC000001", ...)
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

    def save_data(
        self,
        file: str | PathLike,
        links: LinkGraph | None = None,
    ) -> None:
        """Pickle data to a file.

        The  pickled data is a tuple of BGCs, GCFs, Spectra, MolecularFamilies, StrainCollection and
        links, i.e. `(bgcs, gcfs, spectra, mfs, strains, links)`.

        Args:
            file: The path to the pickle file to save the data to.
            links: The LinkGraph object to save.

        Examples:
            Saving the data to a pickle file, links data is `None`:
            >>> npl.save_data("path/to/output.pkl")

            Also saving the links data:
            >>> lg = npl.get_links(npl.gcfs, "metcalf")
            >>> npl.save_data("path/to/output.pkl", lg)
        """
        data = (self.bgcs, self.gcfs, self.spectra, self.mfs, self.strains, links)
        with open(file, "wb") as f:
            pickle.dump(data, f)

    def export_objects(self, objects: Sequence[BGC | Spectrum], filename: str) -> None:
        """Exports the data for a list of BGC or Spectrum objects to a specified file in tab-separated format.

        Args:
            objects (list[BGC | Spectrum]): A list of BGC or Spectrum objects to be exported.
            filename (str): The name of the file where the data will be saved.
        """
        headers = objects[0].to_dict().keys()
        with open(self._output_dir / filename, "w") as f:
            f.write("\t".join(headers) + "\n")
            for obj in objects:
                row_data = obj.to_dict()
                formatted_row = []
                for header in headers:
                    item = row_data.get(header, "")
                    # Convert list, tuple, set to comma-separated string
                    if isinstance(item, (list, tuple, set)):
                        formatted_row.append(", ".join(map(str, item)))
                    # Convert dict to comma-separated string
                    elif isinstance(item, dict):
                        formatted_row.append(", ".join([f"{k}:{v}" for k, v in item.items()]))
                    # Convert non-empty value to string
                    elif item:
                        formatted_row.append(str(item))
                    # Convert empty value to empty string
                    else:
                        formatted_row.append("")
                f.write("\t".join(formatted_row) + "\n")

    def export_results(self, lg: LinkGraph | None = None) -> None:
        """Exports the results to the output directory in tab-separated format.

        This method exports genomics and metabolomics data to their respective
        TSV files in the specified output directory. If a LinkGraph object is
        provided, it also exports the links data to a TSV file.

        Args:
            lg (LinkGraph | None): An optional LinkGraph object. If provided,
                       the links data will be exported to 'links.tsv'.
        """
        self.export_objects(self.bgcs, "genomics_data.tsv")
        self.export_objects(self.spectra, "metabolomics_data.tsv")
        if lg is not None:
            lg.export_links(self._output_dir / "links.tsv")
