from __future__ import annotations
import csv
from os import PathLike
from nplinker.metabolomics.abc import MolecularFamilyLoaderBase
from nplinker.utils import is_file_format
from ..molecular_family import MolecularFamily


class GNPSMolecularFamilyLoader(MolecularFamilyLoaderBase):
    """Load molecular families from GNPS data.

    ??? info "Concept"
        [GNPS data][gnps-data]

    The molecular family file is from GNPS output archive, as described below
    for each GNPS workflow type:

    1. METABOLOMICS-SNETS
        - networkedges_selfloop/*.pairsinfo
    2. METABOLOMICS-SNETS-V2
        - networkedges_selfloop/*.selfloop
    3. FEATURE-BASED-MOLECULAR-NETWORKING
        - networkedges_selfloop/*.selfloop

    The `ComponentIndex` column in the GNPS molecular family file is treated
    as family id.

    But for molecular families that have only one member (i.e. spectrum),
    named singleton molecular families, their files have the same value of
    `-1` in the `ComponentIndex` column. To make the family id unique,the
    spectrum id plus a prefix `singleton-` is used as the family id of
    singleton molecular families.
    """

    def __init__(self, file: str | PathLike) -> None:
        """Initialize the GNPSMolecularFamilyLoader.

        Args:
            file: Path to the GNPS molecular family file.

        Raises:
            ValueError: Raises ValueError if the file is not valid.

        Examples:
            >>> loader = GNPSMolecularFamilyLoader("gnps_molecular_families.tsv")
            >>> print(loader.families)
            [<MolecularFamily 1>, <MolecularFamily 2>, ...]
            >>> print(loader.families[0].spectra_ids)
            {'1', '3', '7', ...}
        """
        self._mfs: list[MolecularFamily] = []
        self._file = file

        self._validate()
        self._load()

    def get_mfs(self, keep_singleton: bool = False) -> list[MolecularFamily]:
        """Get MolecularFamily objects.

        Args:
            keep_singleton: True to keep singleton molecular families. A
                singleton molecular family is a molecular family that contains
                only one spectrum.

        Returns:
            A list of MolecularFamily objects with their spectra ids.
        """
        mfs = self._mfs
        if not keep_singleton:
            mfs = [mf for mf in mfs if not mf.is_singleton()]
        return mfs

    def _validate(self) -> None:
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
                mf_id = row["ComponentIndex"]
                if mf_id not in family_dict:
                    family_dict[mf_id] = set([spec1_id, spec2_id])
                else:
                    family_dict[mf_id].add(spec1_id)
                    family_dict[mf_id].add(spec2_id)
        # convert dict to list of MolecularFamily objects
        for mf_id, spectra_ids in family_dict.items():
            if mf_id == "-1":  # "-1" is from GNPS, it means the singleton molecular family
                for spectrum_id in spectra_ids:
                    # family id must be unique, so using "singleton-" + spectrum id as family id
                    family = MolecularFamily("singleton-" + str(spectrum_id))
                    family.spectra_ids = set([spectrum_id])
                    self._mfs.append(family)
            else:
                # for regular molecular families, use the value of "ComponentIndex" as family id
                family = MolecularFamily(mf_id)
                family.spectra_ids = spectra_ids
                self._mfs.append(family)
