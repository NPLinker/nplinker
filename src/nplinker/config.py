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
    Validator("loglevel", is_type_of=str),
    Validator("logfile", is_type_of=str),
    Validator("repro_file", is_type_of=str),
    Validator("log_to_stdout", is_type_of=bool),
]
config.validators.register(*validators)
config.validators.validate()
