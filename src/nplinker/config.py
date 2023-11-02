import argparse
import os
from collections.abc import Mapping
from shutil import copyfile
import toml
from xdg import XDG_CONFIG_HOME
from .logconfig import LogConfig


try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

logger = LogConfig.getLogger(__name__)


class Args:
    def __init__(self):
        def bool_checker(x):
            return str(x).lower() == "true"

        # TODO need to finalise these
        self.parser = argparse.ArgumentParser(
            description="nplinker arguments",
            epilog="Note: command-line arguments will override "
            "arguments from configuration files",
        )
        self.parser.add_argument(
            "-c", "--config", help="Path to a .toml configuration file", metavar="path"
        )
        self.parser.add_argument(
            "-d",
            "--dataset.root",
            help='Root path for the dataset to be loaded (or "platform:datasetID" for remote datasets)',
            metavar="root",
        )
        self.parser.add_argument(
            "-l",
            "--loglevel",
            help="Logging verbosity level: DEBUG, INFO, WARNING, ERROR",
            metavar="loglevel",
        )
        self.parser.add_argument(
            "-f", "--logfile", help="Redirect logging from stdout to this file", metavar="logfile"
        )
        self.parser.add_argument(
            "-s",
            "--log_to_stdout",
            help="keep logging to stdout even if --logfile used",
            metavar="log_to_stdout",
        )

        self.parser.add_argument(
            "--bigscape-cutoff", help="BIGSCAPE clustering cutoff threshold", metavar="cutoff"
        )
        self.parser.add_argument(
            "--repro-file", help="Filename to store reproducibility data in", metavar="filename"
        )

        self.args = self.parser.parse_args()

    def get_args(self):
        # restructure the arguments into a dict with the same nested structure as Config class expects
        orig = vars(self.args)
        args = {}
        for k, v in orig.items():
            # ignore any params with no value given
            if v is None:
                continue

            # values with non-dotted names can get inserted directly
            if k.find(".") == -1:
                args[k] = v
            else:
                # otherwise add a nested dict for each dotted part, then
                # insert the actual value on the innermost level
                parts = k.split(".")
                root = args
                for p in parts[:-1]:
                    if p not in root:
                        root[p] = {}
                    root = root[p]
                root[parts[-1]] = v
        return args


class Config:
    """Wrapper for all NPLinker configuration options"""

    DEFAULT_CONFIG = "nplinker.toml"

    def __init__(self, config_dict):
        self.default_config_path = os.path.join(XDG_CONFIG_HOME, "nplinker", Config.DEFAULT_CONFIG)
        if not os.path.exists(self.default_config_path):
            logger.debug("Creating default config file")
            os.makedirs(os.path.join(XDG_CONFIG_HOME, "nplinker"), exist_ok=True)
            copyfile(
                files("nplinker").joinpath("data", Config.DEFAULT_CONFIG), self.default_config_path
            )

        # load the default per-user config file, then check for one provided as an argument
        # and if present use it to override the defaults
        logger.debug("Parsing default config file: {}".format(self.default_config_path))
        config = toml.load(open(self.default_config_path))
        if "config" in config_dict:
            logger.debug("Loading user config {}".format(config_dict["config"]))
            if config_dict["config"] is not None:
                user_config = toml.load(open(config_dict["config"]))
                config.update(user_config)
                del config_dict["config"]

        # remaining values in the dict should override the existing ones from config files
        # however if running non-interactively, argparse will set values of all non-specified
        # options to None and don't want to wipe out existing settings, so do things this way
        def update(d, u):
            for k, v in u.items():
                if isinstance(v, Mapping):
                    d[k] = update(d.get(k, {}), v)
                elif v is not None:
                    d[k] = v
            return d

        config = update(config, config_dict)
        self._validate(config)
        self.config = config

    def _validate(self, config: dict) -> None:
        """Validates the configuration dictionary to ensure that all required
            fields are present and have valid values.

        Args:
            config (dict): The configuration dictionary to validate.

        Raises:
            ValueError: If the configuration dictionary is missing required
                fields or contains invalid values.
        """
        if "dataset" not in config:
            raise ValueError('Not found config for "dataset".')

        root = config["dataset"].get("root")
        if root is None:
            raise ValueError('Not found config for "root".')

        if root.startswith("platform:"):
            config["dataset"]["platform_id"] = root.replace("platform:", "")
            logger.info("Loading from platform project ID %s", config["dataset"]["platform_id"])
        else:
            config["dataset"]["platform_id"] = ""
            logger.info("Loading from local data in directory %s", root)

        antismash = config["dataset"].get("antismash")
        allowed_antismash_formats = ["default", "flat"]
        if antismash is not None:
            if "format" in antismash and antismash["format"] not in allowed_antismash_formats:
                raise ValueError(f'Unknown antismash format: {antismash["format"]}')
