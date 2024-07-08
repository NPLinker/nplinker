from __future__ import annotations
import csv
from os import PathLike
from pathlib import Path
from nplinker.metabolomics.abc import FileMappingLoaderBase
from nplinker.utils import is_file_format
from .gnps_format import GNPSFormat
from .gnps_format import gnps_format_from_file_mapping


class GNPSFileMappingLoader(FileMappingLoaderBase):
    """Class to load file mappings from GNPS output file.

    ??? info "Concept"
        [GNPS data][gnps-data]

    File mappings refers to the mapping from spectrum id to files in which
    this spectrum occurs.

    The file mappings file is from GNPS output archive, as described below
    for each GNPS workflow type:

    1. METABOLOMICS-SNETS
        - clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.tsv
    2. METABOLOMICS-SNETS-V2
        - clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.clustersummary
    3. FEATURE-BASED-MOLECULAR-NETWORKING
        - quantification_table*/*.csv
    """

    def __init__(self, file: str | PathLike) -> None:
        """Initialize the GNPSFileMappingLoader.

        Args:
            file: Path to the GNPS file mappings file.

        Raises:
            ValueError: Raises ValueError if the file is not valid.

        Examples:
            >>> loader = GNPSFileMappingLoader("gnps_file_mappings.tsv")
            >>> print(loader.mappings["1"])
            ['26c.mzXML']
            >>> print(loader.mapping_reversed["26c.mzXML"])
            {'1', '3', '7', ...}
        """
        self._gnps_format = gnps_format_from_file_mapping(file)
        if self._gnps_format is GNPSFormat.Unknown:
            raise ValueError("Unknown workflow type for GNPS file mappings file ")

        self._file = Path(file)
        self._mapping: dict[str, list[str]] = {}

        self._validate()
        self._load()

    @property
    def mappings(self) -> dict[str, list[str]]:
        """Return mapping from spectrum id to files in which this spectrum occurs.

        Returns:
            Mapping from spectrum id to names of all files in which this spectrum occurs.
        """
        return self._mapping

    @property
    def mapping_reversed(self) -> dict[str, set[str]]:
        """Return mapping from file name to all spectra that occur in this file.

        Returns:
            Mapping from file name to all spectra ids that occur in this file.
        """
        mapping_reversed: dict[str, set[str]] = {}
        for spectrum_id, ms_filenames in self._mapping.items():
            for filename in ms_filenames:
                if filename in mapping_reversed:
                    mapping_reversed[filename].add(spectrum_id)
                else:
                    mapping_reversed[filename] = {spectrum_id}

        return mapping_reversed

    def _validate(self) -> None:
        """Validate the file mappings file.

        Raises:
            ValueError: Raises ValueError if the file is not valid.
        """
        # validate file format
        required_file_formats = {
            GNPSFormat.SNETS: "tsv",
            GNPSFormat.SNETSV2: "tsv",
            GNPSFormat.FBMN: "csv",
        }
        if not is_file_format(self._file, required_file_formats[self._gnps_format]):
            raise ValueError(
                f"Invalid GNPS file mappings file '{self._file}'. "
                f"Expected a {required_file_formats[self._gnps_format]} file."
            )

        # validate required columns against the header
        required_columns = {
            GNPSFormat.SNETS: ["cluster index", "AllFiles"],
            GNPSFormat.SNETSV2: ["cluster index", "UniqueFileSources"],
            GNPSFormat.FBMN: ["row ID", " Peak area"],
        }
        with open(self._file, mode="rt") as f:
            header = f.readline()
            for k in required_columns[self._gnps_format]:
                if k not in header:
                    raise ValueError(
                        f"Invalid GNPS file mappings file '{self._file}'. "
                        f"Expected a header line with '{k}' column, "
                        f"but got '{header}'."
                    )

        # validate that cluster index or row id must be unique
        with open(self._file, mode="rt") as f:
            if self._gnps_format is GNPSFormat.FBMN:
                reader = csv.DictReader(f, delimiter=",")
                ids = [row["row ID"] for row in reader]
            else:
                reader = csv.DictReader(f, delimiter="\t")
                ids = [row["cluster index"] for row in reader]
        duplicates = {x for x in ids if ids.count(x) > 1}
        if len(duplicates) > 0:
            raise ValueError(
                f"Invalid GNPS file mappings file '{self._file}'. "
                f"Expected unique 'cluster index' or 'row ID', "
                f"but found duplicates '{duplicates}'."
            )

    def _load(self) -> None:
        """Load file mapping from the file based on the GNPS workflow type."""
        if self._gnps_format is GNPSFormat.SNETS:
            self._load_snets()
        elif self._gnps_format is GNPSFormat.SNETSV2:
            self._load_snetsv2()
        elif self._gnps_format is GNPSFormat.FBMN:
            self._load_fbmn()

    def _load_snets(self) -> None:
        """Load file mapping from output of GNPS SNETS workflow.

        The following columns are loaded:

        - "cluster index": loaded as spectrum id
        - "AllFiles": a list of files in which the spectrum occurs, separated
            by '###'.
            An example data of "AllFiles" column is as follows:
            "2b.mzXML:1503195###6a.mzXML:1502983###"
        """
        with open(self._file, mode="rt", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                spectrum_id = row["cluster index"]
                occurrences = row["AllFiles"].split("###")  # split by '###'
                occurrences.pop()  # remove last empty entry
                # separate the scan position from the files
                samples = [x.split(":")[0] for x in occurrences]
                self._mapping[spectrum_id] = samples

    def _load_snetsv2(self) -> None:
        """Load file mapping from output of GNPS SNETS-V2 workflow.

        The following columns are loaded:

        - "cluster index": loaded as spectrum id
        - "UniqueFileSources": a list of files in which the spectrum occurs,
            separated by '|'.
            An example data of "UniqueFileSources" column is as follows:
            "140221_Blanc5.mzML|140221_Blanc8.mzML|140221_ME_14_12.mzML"
        """
        with open(self._file, mode="rt", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                spectrum_id = row["cluster index"]
                samples = row["UniqueFileSources"].split("|")
                self._mapping[spectrum_id] = samples

    def _load_fbmn(self) -> None:
        """Load file mapping from output of GNPS FBMN workflow.

        The column "row ID" is loaded as spectrum id.

        The column names containing " Peak area" are used to extract the file
        names, and the values of these columns are used to determine whether
        the spectrum occurs in the file. The file name is taken only if the
        value is greater than 0.

        An example data of the file is as follows:
            ```
            row ID,5434_5433_mod.mzXML Peak area,5425_5426_mod.mzXML Peak area
            1,1764067.8434999974,0.0
            ```
        """
        pattern = " Peak area"
        with open(self._file, mode="rt", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=",")
            for row in reader:
                spectrum_id = row["row ID"]
                samples = []
                for col in row:
                    if pattern in col and float(row[col]) > 0:
                        samples.append(col.strip(pattern))
                self._mapping[spectrum_id] = samples
