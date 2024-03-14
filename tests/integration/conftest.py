import os
import tempfile
import zipfile
from . import DATA_DIR


# Prepare dataset for local mode testing
# ⚠️ Multiple temp dirs will be created if using parallel testing.
temp_dir = tempfile.mkdtemp(prefix="nplinker_")
nplinker_root_dir = os.path.join(temp_dir, "local_mode_example")
with zipfile.ZipFile(DATA_DIR / "local_mode_example.zip", "r") as zip_ref:
    zip_ref.extractall(temp_dir)

# NPLinker setting `root_dir` must be a path that exists, so setting it to a temporary directory.
os.environ["NPLINKER_ROOT_DIR"] = nplinker_root_dir

# Specify the config file via environment variable before importing nplinker in any test.
os.environ["NPLINKER_CONFIG_FILE"] = str(DATA_DIR / "nplinker_local_mode.toml")
