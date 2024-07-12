import json
from pathlib import Path
from jsonschema import Draft7Validator


__all__ = [
    "GENOME_STATUS_SCHEMA",
    "GENOME_BGC_MAPPINGS_SCHEMA",
    "STRAIN_MAPPINGS_SCHEMA",
    "PODP_ADAPTED_SCHEMA",
    "USER_STRAINS_SCHEMA",
    "validate_podp_json",
]

_schema_dir = Path(__file__).parent

with open(_schema_dir / "genome_status_schema.json", "r") as f:
    GENOME_STATUS_SCHEMA = json.load(f)
    """Schema for the [genome status JSON file][nplinker.defaults.GENOME_STATUS_FILENAME].

    ??? info "Schema Content: `genome_status_schema.json`"
        ```json
        --8<-- "src/nplinker/schemas/genome_status_schema.json"
        ```
    """

with open(_schema_dir / "genome_bgc_mappings_schema.json", "r") as f:
    GENOME_BGC_MAPPINGS_SCHEMA = json.load(f)
    """Schema for [genome BGC mappings JSON file][nplinker.defaults.GENOME_BGC_MAPPINGS_FILENAME].

    ??? info "Schema Content: `genome_bgc_mappings_schema.json`"
        ```json
        --8<-- "src/nplinker/schemas/genome_bgc_mappings_schema.json"
        ```
    """

with open(_schema_dir / "strain_mappings_schema.json", "r") as f:
    STRAIN_MAPPINGS_SCHEMA = json.load(f)
    """Schema for [strain mappings JSON file][nplinker.defaults.STRAIN_MAPPINGS_FILENAME].

    ??? info "Schema Content: `strain_mappings_schema.json`"
        ```json
        --8<-- "src/nplinker/schemas/strain_mappings_schema.json"
        ```
    """

with open(_schema_dir / "user_strains.json", "r") as f:
    USER_STRAINS_SCHEMA = json.load(f)
    """Schema for [user strains JSON file][nplinker.defaults.STRAINS_SELECTED_FILENAME].

    ??? info "Schema Content: `user_strains.json`"
        ```json
        --8<-- "src/nplinker/schemas/user_strains.json"
        ```
    """

with open(_schema_dir / "podp_adapted_schema.json", "r") as f:
    PODP_ADAPTED_SCHEMA = json.load(f)
    """Schema for PODP JSON file.

    The PODP JSON file is the project JSON file downloaded from PODP platform.
    For example, for PODP project MSV000079284, its JSON file is
    https://pairedomicsdata.bioinformatics.nl/api/projects/4b29ddc3-26d0-40d7-80c5-44fb6631dbf9.4.

    ??? info "Schema Content: `podp_adapted_schema.json`"
        ```json
        --8<-- "src/nplinker/schemas/podp_adapted_schema.json"
        ```
    """


def validate_podp_json(json_data: dict) -> None:
    """Validate JSON data against the PODP JSON schema.

    All validation error messages are collected and raised as a single
    ValueError.

    Args:
        json_data: The JSON data to validate.

    Raises:
        ValueError: If the JSON data does not match the schema.

    Examples:
        Download PODP JSON file for project MSV000079284 from
        https://pairedomicsdata.bioinformatics.nl/api/projects/4b29ddc3-26d0-40d7-80c5-44fb6631dbf9.4
        and save it as `podp_project.json`.

        Validate it:
        >>> with open(podp_project.json, "r") as f:
        ...     json_data = json.load(f)
        >>> validate_podp_json(json_data)
    """
    validator = Draft7Validator(PODP_ADAPTED_SCHEMA)
    errors = sorted(validator.iter_errors(json_data), key=lambda e: e.path)
    if errors:
        error_messages = [f"{e.json_path}: {e.message}" for e in errors]
        raise ValueError(
            "Not match PODP adapted schema, here are the detailed error:\n  - "
            + "\n  - ".join(error_messages)
        )
