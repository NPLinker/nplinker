import os
from . import DATA_DIR


# Setting the environment variable to the config file path before importing nplinker in any test
os.environ["NPLINKER_CONFIG_FILE"] = str(DATA_DIR / "nplinker_demo1.toml")
