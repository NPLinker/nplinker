import os
import tempfile
from . import DATA_DIR


# Prepare root dir
nplinker_root_dir = os.path.join(tempfile.gettempdir(), "nplinker_unit_test")
if not os.path.exists(nplinker_root_dir):
    os.mkdir(nplinker_root_dir)

# NPLinker setting `root_dir` must be a path that exists, so setting it to a temporary directory.
os.environ["NPLINKER_ROOT_DIR"] = nplinker_root_dir

# Specify the config file via environment variable before importing nplinker in any test.
os.environ["NPLINKER_CONFIG_FILE"] = str(DATA_DIR / "nplinker_local_mode.toml")
