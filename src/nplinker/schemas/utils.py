import json
from pathlib import Path
from jsonschema import Draft7Validator


with open(Path(__file__).parent / "podp_adapted_schema.json", 'r') as f:
    PODP_ADAPTED_SCHEMA = json.load(f)


def validate_podp_json(json_data: dict) -> None:
    """
    Validate a dictionary of JSON data against the PODP JSON schema.

    All validation error messages are collected and raised as a single
    ValueError.

    Parameters:
        json_data (dict): The JSON data to validate.

    Raises:
        ValueError: If the JSON data does not match the schema.
    """
    validator = Draft7Validator(PODP_ADAPTED_SCHEMA)
    errors = sorted(validator.iter_errors(json_data), key=lambda e: e.path)
    if errors:
        error_messages = [f"{e.json_path}: {e.message}" for e in errors]
        raise ValueError(
            "Not match PODP adapted schema, here are the detailed error:\n  - "
            + "\n  - ".join(error_messages))
