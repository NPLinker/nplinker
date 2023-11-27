import csv
from os import PathLike
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import SingletonFamily
from nplinker.metabolomics.abc import MolecularFamilyLoaderBase
from nplinker.utils import is_file_format


class GNPSMolecularFamilyLoader(MolecularFamilyLoaderBase):
    def __init__(self, file: str | PathLike):
        """Class to load molecular families from GNPS output file.

        The molecular family file is from GNPS output archive, as described below
        for each GNPS workflow type:
        1. METABOLOMICS-SNETS
            - networkedges_selfloop/*.pairsinfo
        2. METABOLOMICS-SNETS-V2
            - networkedges_selfloop/*.selfloop
        3. FEATURE-BASED-MOLECULAR-NETWORKING
            - networkedges_selfloop/*.selfloop

        Args:
            file(str | PathLike): Path to the GNPS molecular family file.

        Raises:
            ValueError: Raises ValueError if the file is not valid.

        Example:
            >>> loader = GNPSMolecularFamilyLoader("gnps_molecular_families.tsv")
            >>> print(loader.families)
            [<MolecularFamily 1>, <MolecularFamily 2>, ...]
            >>> print(loader.families[0].spectra_ids)
            {'1', '3', '7', ...}
        """
        self._mfs: list[MolecularFamily | SingletonFamily] = []
        self._file = file

        self._validate()
        self._load()

    def get_mfs(self, keep_singleton: bool = False) -> list[MolecularFamily]:
        """Get MolecularFamily objects.

        Args:
            keep_singleton(bool): True to keep singleton molecular families. A
                singleton molecular family is a molecular family that contains
                only one spectrum.

        Returns:
            list[MolecularFamily]: A list of MolecularFamily objects with their
                spectra ids.
        """
        mfs = self._mfs
        if not keep_singleton:
            mfs = [mf for mf in mfs if not mf.is_singleton()]
        return mfs

    def _validate(self):
        """Validate the GNPS molecular family file."""
        # validate file format
        if not is_file_format(self._file, "tsv"):
            raise ValueError(
                f"Invalid GNPS molecular family file '{self._file}'. " f"Expected a '.tsv' file."
            )
        # validate required columns against the header
        required_columns = ["CLUSTERID1", "CLUSTERID2", "ComponentIndex"]
        with open(self._file, mode="rt") as f:
            header = f.readline()
            for k in required_columns:
                if k not in header:
                    raise ValueError(
                        f"Invalid GNPS molecular famliy file '{self._file}'. "
                        f"Expected a header line with '{k}' column, "
                        f"but got '{header}'."
                    )

    def _load(self) -> None:
        """Load molecular families from GNPS output file.

        Molecular families are loaded as a list of MolecularFamily objects. Each
        MolecularFamily object contains a set of spectra ids that belong to this
        family.
        """
        # load molecular families to dict
        family_dict = {}
        with open(self._file, mode="rt", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                spec1_id = row["CLUSTERID1"]
                spec2_id = row["CLUSTERID2"]
                family_id = row["ComponentIndex"]
                if family_id not in family_dict:
                    family_dict[family_id] = set([spec1_id, spec2_id])
                else:
                    family_dict[family_id].add(spec1_id)
                    family_dict[family_id].add(spec2_id)
        # convert dict to list of MolecularFamily objects
        for family_id, spectra_ids in family_dict.items():
            if family_id == "-1":  # the "-1" is from GNPS result
                for spectrum_id in spectra_ids:
                    family = SingletonFamily()  ## uuid as family id
                    family.spectra_ids = set([spectrum_id])
                    self._mfs.append(family)
            else:
                family = MolecularFamily(family_id)
                family.spectra_ids = spectra_ids
                self._mfs.append(family)
