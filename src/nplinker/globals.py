from pathlib import Path
from nplinker.config import config


STRAIN_MAPPINGS_FILENAME = "strain_mappings.json"
GENOME_BGC_MAPPINGS_FILENAME = "genome_bgc_mappings.json"
GENOME_STATUS_FILENAME = "genome_status.json"
GNPS_FILE_MAPPINGS_FILENAME = "file_mappings"


DOWNLOADS_DEFAULT_PATH: Path = config.root_dir / "downloads"
MIBIG_DEFAULT_PATH: Path = config.root_dir / "mibig"
GNPS_DEFAULT_PATH: Path = config.root_dir / "gnps"
ANTISMASH_DEFAULT_PATH: Path = config.root_dir / "antismash"
BIGSCAPE_DEFAULT_PATH: Path = config.root_dir / "bigscape"
