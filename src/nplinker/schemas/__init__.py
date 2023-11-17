import json
import logging
from pathlib import Path
from .utils import PODP_ADAPTED_SCHEMA
from .utils import validate_podp_json


logging.getLogger(__name__).addHandler(logging.NullHandler())

__all__ = [
    "GENOME_STATUS_SCHEMA",
    "GENOME_BGC_MAPPINGS_SCHEMA",
    "STRAIN_MAPPINGS_SCHEMA",
    "PODP_ADAPTED_SCHEMA",
    "USER_STRAINS_SCHEMA",
    "validate_podp_json",
]

SCHEMA_DIR = Path(__file__).parent
with open(SCHEMA_DIR / "genome_status_schema.json", "r") as f:
    GENOME_STATUS_SCHEMA = json.load(f)

with open(SCHEMA_DIR / "genome_bgc_mappings_schema.json", "r") as f:
    GENOME_BGC_MAPPINGS_SCHEMA = json.load(f)

with open(SCHEMA_DIR / "strain_mappings_schema.json", "r") as f:
    STRAIN_MAPPINGS_SCHEMA = json.load(f)

with open(SCHEMA_DIR / "user_strains.json", "r") as f:
    USER_STRAINS_SCHEMA = json.load(f)
