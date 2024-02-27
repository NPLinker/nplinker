import os
import tempfile
from . import DATA_DIR


# Specify the config file via environment variable before importing nplinker in any test.
os.environ["NPLINKER_CONFIG_FILE"] = str(DATA_DIR / "nplinker_demo1.toml")
# NPLinker setting `root_dir` must be a path that exists, so setting it to a temporary directory.
os.environ["NPLINKER_ROOT_DIR"] = tempfile.mkdtemp(prefix="nplinker_")
