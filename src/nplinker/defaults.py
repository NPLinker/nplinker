from pathlib import Path
from typing import Final


# The path to the NPLinker application database directory
NPLINKER_APP_DATA_DIR: Final = Path(__file__).parent / "data"


STRAIN_MAPPINGS_FILENAME: Final = "strain_mappings.json"
GENOME_BGC_MAPPINGS_FILENAME: Final = "genome_bgc_mappings.json"
GENOME_STATUS_FILENAME: Final = "genome_status.json"
GNPS_SPECTRA_FILENAME: Final = "spectra.mgf"
GNPS_MOLECULAR_FAMILY_FILENAME: Final = "molecular_families.tsv"
GNPS_ANNOTATIONS_FILENAME: Final = "annotations.tsv"
GNPS_FILE_MAPPINGS_TSV: Final = "file_mappings.tsv"
GNPS_FILE_MAPPINGS_CSV: Final = "file_mappings.csv"
STRAINS_SELECTED_FILENAME: Final = "strains_selected.json"


DOWNLOADS_DIRNAME: Final = "downloads"
MIBIG_DIRNAME: Final = "mibig"
GNPS_DIRNAME: Final = "gnps"
ANTISMASH_DIRNAME: Final = "antismash"
BIGSCAPE_DIRNAME: Final = "bigscape"
BIGSCAPE_RUNNING_OUTPUT_DIRNAME: Final = "bigscape_running_output"
OUTPUT_DIRNAME: Final = "output"
