import os
from pathlib import Path
from dynaconf import Dynaconf
from dynaconf import Validator


# The Dynaconf library is used for loading NPLinker config file.
# Users can set the config file location via the NPLINKER_CONFIG_FILE environment variable before
# running/importing NPLinker. If not set, we default to 'nplinker.toml' file in the current working
# directory.
# The loaded config data is available by importing this module and accessing the 'config' variable.


# Locate the user's config file
user_config_file = os.environ.get("NPLINKER_CONFIG_FILE", "nplinker.toml")
if not os.path.exists(user_config_file):
    raise FileNotFoundError(f"Config file '{user_config_file}' not found")

# Locate the default config file
default_config_file = Path(__file__).resolve().parent / "nplinker_default.toml"

# Load config files
config = Dynaconf(settings_files=[user_config_file], preload=[default_config_file])

# Validate config
# Note:
# Validataor parameter `required=False` means the setting (e.g. "loglevel") must not exist rather
# than being optional. So don't set the parameter `required` if the key is optional.
validators = [
    # General settings
    ## `root_dir` value is transformed to a `pathlib.Path` object and must be a directory.
    Validator("root_dir", required=True, cast=Path, condition=lambda v: v.is_dir()),
    Validator("mode", required=True, cast=lambda v: v.lower(), is_in=["local", "podp"]),
    ## `podp_id` must be set if `mode` is "podp"; must not be set if `mode` is "local".
    Validator("podp_id", required=True, when=Validator("mode", eq="podp")),
    Validator("podp_id", required=False, when=Validator("mode", eq="local")),
    # Log
    ## `loglevel` must be a string and must be one of the supported levels. It is transformed to
    ## uppercase to avoid case sensitivity.
    Validator(
        "log.level",
        is_type_of=str,
        cast=lambda v: v.upper(),
        is_in=["NOTSET", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    ),
    Validator("log.file", is_type_of=str, cast=Path),
    Validator("log.to_stdout", is_type_of=bool),
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
config.validators.register(*validators)
config.validators.validate()
