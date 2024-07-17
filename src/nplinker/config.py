from __future__ import annotations
from os import PathLike
from pathlib import Path
from dynaconf import Dynaconf
from dynaconf import Validator
from nplinker.utils import transform_to_full_path


def load_config(config_file: str | PathLike) -> Dynaconf:
    """Load and validate the configuration file.

    ??? info "Usage Documentation"
        [Config Loader][config-loader]

    Args:
        config_file: Path to the configuration file.

    Returns:
        Dynaconf: A Dynaconf object containing the configuration settings.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
    """
    config_file = transform_to_full_path(config_file)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file '{config_file}' not found")

    # Locate the default config file
    default_config_file = Path(__file__).resolve().parent / "nplinker_default.toml"

    # Load config files
    config = Dynaconf(settings_files=[config_file], preload=[default_config_file])

    # Validate configs
    config.validators.register(*CONFIG_VALIDATORS)
    config.validators.validate()

    return config


# Note:
# Validator parameter `required=False` means the setting (e.g. "loglevel") must not exist rather
# than being optional. So don't set the parameter `required` if the key is optional.
CONFIG_VALIDATORS = [
    # General settings
    ## `root_dir` value is transformed to a `pathlib.Path` object and must be a directory.
    Validator(
        "root_dir", required=True, cast=transform_to_full_path, condition=lambda v: v.is_dir()
    ),
    Validator("mode", required=True, cast=lambda v: v.lower(), is_in=["local", "podp"]),
    ## `podp_id` must be set if `mode` is "podp"; must not be set if `mode` is "local".
    Validator("podp_id", required=True, when=Validator("mode", eq="podp")),
    Validator("podp_id", required=False, when=Validator("mode", eq="local")),
    # Log
    ## `level` must be a string and must be one of the supported levels. It is transformed to
    ## uppercase to avoid case sensitivity.
    Validator(
        "log.level",
        is_type_of=str,
        cast=lambda v: v.upper(),
        is_in=["NOTSET", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    ),
    Validator("log.file", is_type_of=str),
    Validator("log.use_console", is_type_of=bool),
    #  Mibig
    Validator("mibig.to_use", required=True, is_type_of=bool),
    Validator(
        "mibig.version",
        required=True,
        is_type_of=str,
        when=Validator("mibig.to_use", eq=True),
    ),
    # BigScape
    Validator("bigscape.parameters", required=True, is_type_of=str),
    Validator("bigscape.cutoff", required=True, is_type_of=str),
    Validator("bigscape.version", required=True, is_type_of=int),
    # Scoring
    ## `scoring.methods` must be a list of strings and must contain at least one of the
    ## supported scoring methods.
    Validator(
        "scoring.methods",
        required=True,
        cast=lambda v: [i.lower() for i in v],
        is_type_of=list,
        len_min=1,
        condition=lambda v: set(v).issubset({"metcalf", "rosetta"}),
    ),
]
