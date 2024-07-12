from __future__ import annotations
import csv
from os import PathLike
from pathlib import Path
from nplinker.metabolomics.abc import AnnotationLoaderBase
from nplinker.utils import is_file_format


GNPS_UNIVERSAL_SPECTRUM_IDENTIFIER_URL = (
    "https://metabolomics-usi.gnps2.org/{}/?usi1=mzspec:GNPS:GNPS-LIBRARY:accession:{}"
)


class GNPSAnnotationLoader(AnnotationLoaderBase):
    """Load annotations from GNPS output file.

    ??? info "Concept"
        [GNPS data][gnps-data]

    The annotation file is a `.tsv` file from GNPS output archive, as described
    below for each GNPS workflow type:

    1. METABOLOMICS-SNETS
        - result_specnets_DB/*.tsv
    2. METABOLOMICS-SNETS-V2
        - result_specnets_DB/.tsv
    3. FEATURE-BASED-MOLECULAR-NETWORKING
        - DB_result/*.tsv
    """

    def __init__(self, file: str | PathLike) -> None:
        """Initialize the GNPSAnnotationLoader.

        Args:
            file: The GNPS annotation file.

        Examples:
            >>> loader = GNPSAnnotationLoader("gnps_annotations.tsv")
            >>> print(loader.annotations["100"])
            {'#Scan#': '100',
            'Adduct': 'M+H',
            'CAS_Number': 'N/A',
            'Charge': '1',
            'Compound_Name': 'MLS002153841-01!Iobenguane sulfate',
            'Compound_Source': 'NIH Pharmacologically Active Library',
            'Data_Collector': 'VP/LMS',
            'ExactMass': '274.992',
            'INCHI': 'N/A',
            'INCHI_AUX': 'N/A',
            'Instrument': 'qTof',
            'IonMode': 'Positive',
            'Ion_Source': 'LC-ESI',
            'LibMZ': '276.003',
            'LibraryName': 'lib-00014.mgf',
            'LibraryQualityString': 'Gold',
            'Library_Class': '1',
            'MQScore': '0.704152',
            'MZErrorPPM': '405416',
            'MassDiff': '111.896',
            'Organism': 'GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE',
            'PI': 'Dorrestein',
            'Precursor_MZ': '276.003',
            'Pubmed_ID': 'N/A',
            'RT_Query': '795.979',
            'SharedPeaks': '7',
            'Smiles': 'NC(=N)NCc1cccc(I)c1.OS(=O)(=O)O',
            'SpecCharge': '1',
            'SpecMZ': '164.107',
            'SpectrumFile': 'spectra/specs_ms.pklbin',
            'SpectrumID': 'CCMSLIB00000086167',
            'TIC_Query': '986.997',
            'UpdateWorkflowName': 'UPDATE-SINGLE-ANNOTATED-GOLD',
            'tags': ' ',
            'png_url': 'https://metabolomics-usi.gnps2.org/png/?usi1=mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000086167',
            'json_url': 'https://metabolomics-usi.gnps2.org/json/?usi1=mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000086167',
            'svg_url': 'https://metabolomics-usi.gnps2.org/svg/?usi1=mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000086167',
            'spectrum_url': 'https://metabolomics-usi.gnps2.org/spectrum/?usi1=mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000086167'}
        """
        self._file = Path(file)
        self._annotations: dict[str, dict] = {}

        self._validate()
        self._load()

    @property
    def annotations(self) -> dict[str, dict]:
        """Get annotations.

        Returns:
            Keys are spectrum ids ("#Scan#" in annotation file) and values are the annotations dict
            for each spectrum.
        """
        return self._annotations

    def _validate(self) -> None:
        """Validate the annotation file.

        Raises:
            ValueError: Raises ValueError if the file is not valid.
        """
        # validate file format
        if not is_file_format(self._file, "tsv"):
            raise ValueError(
                f"Invalid GNPS annotation file '{self._file}'. " f"Expected a .tsv file."
            )

        # validate required columns against the header
        required_columns = ["#Scan#", "Compound_Name", "Organism", "MQScore", "SpectrumID"]
        with open(self._file, mode="rt") as f:
            header = f.readline()
            for k in required_columns:
                if k not in header:
                    raise ValueError(
                        f"Invalid GNPS annotation file '{self._file}'. "
                        f"Expected a header line with '{k}' column, "
                        f"but got '{header}'."
                    )

        # validate that "#Scan#" must be unique
        with open(self._file, mode="rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            scans = [row["#Scan#"] for row in reader]
        duplicates = {x for x in scans if scans.count(x) > 1}
        if len(duplicates) > 0:
            raise ValueError(
                f"Invalid GNPS annotation file '{self._file}'. "
                f"Expected unique '#Scan#', but found duplicates '{duplicates}'."
            )

    def _load(self) -> None:
        """Load the annotations from the file."""
        with open(self._file, mode="rt") as f:
            dict_reader = csv.DictReader(f, delimiter="\t")
            for row in dict_reader:
                scan_id = row["#Scan#"]
                self._annotations[scan_id] = row
                # insert useful URLs
                for t in ["png", "json", "svg", "spectrum"]:
                    self._annotations[scan_id][f"{t}_url"] = (
                        GNPS_UNIVERSAL_SPECTRUM_IDENTIFIER_URL.format(t, row["SpectrumID"])
                    )
