import os
import sys
from pathlib import Path
from nplinker.config import config


PFAM_PATH = os.path.join(sys.prefix, "nplinker_lib")

STRAIN_MAPPINGS_FILENAME = "strain_mappings.json"
GENOME_BGC_MAPPINGS_FILENAME = "genome_bgc_mappings.json"
GENOME_STATUS_FILENAME = "genome_status.json"
GNPS_SPECTRA_FILENAME = "spectra.mgf"
GNPS_MOLECULAR_FAMILY_FILENAME = "molecular_families.tsv"
GNPS_ANNOTATIONS_FILENAME = "annotations.tsv"
GNPS_FILE_MAPPINGS_TSV = "file_mappings.tsv"
GNPS_FILE_MAPPINGS_CSV = "file_mappings.csv"
STRAINS_SELECTED_FILENAME = "strains_selected.json"


DOWNLOADS_DEFAULT_PATH: Path = config.root_dir / "downloads"
MIBIG_DEFAULT_PATH: Path = config.root_dir / "mibig"
GNPS_DEFAULT_PATH: Path = config.root_dir / "gnps"
ANTISMASH_DEFAULT_PATH: Path = config.root_dir / "antismash"
BIGSCAPE_DEFAULT_PATH: Path = config.root_dir / "bigscape"
BIGSCAPE_RUNNING_OUTPUT_PATH: Path = BIGSCAPE_DEFAULT_PATH / "bigscape_running_output"
