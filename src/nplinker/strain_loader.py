import json
from os import PathLike
from jsonschema import validate
from nplinker.logconfig import LogConfig
from nplinker.schemas import USER_STRAINS_SCHEMA
from .strain import Strain


logger = LogConfig.getLogger(__name__)


def load_user_strains(json_file: str | PathLike) -> set[Strain]:
    """Load user specified strains from a JSON file.

    The JSON file must follow the schema defined in "nplinker/schemas/user_strains_schema.json".
    An example content of the JSON file:
        {"strain_ids": ["strain1", "strain2"]}

    Args:
        json_file(str | PathLike): Path to the JSON file containing user specified strains.

    Returns:
        set[Strain]: A set of user specified strains.
    """
    with open(json_file, "r") as f:
        json_data = json.load(f)

    # validate json data
    validate(instance=json_data, schema=USER_STRAINS_SCHEMA)

    strains = set()
    for strain_id in json_data["strain_ids"]:
        strains.add(Strain(strain_id))

    return strains
