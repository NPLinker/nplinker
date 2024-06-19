import os
from nplinker.config import load_config
from . import CONFIG_FILE_LOCAL_MODE


def test_config(tmp_path):
    """Test loading the default config file."""
    os.environ["NPLINKER_ROOT_DIR"] = str(tmp_path)  # Create a temporary root dir for NPLinker
    config = load_config(CONFIG_FILE_LOCAL_MODE)

    assert config.mode == "local"
    assert config.log.level == "DEBUG"
    assert config["log.level"] == "DEBUG"
    assert config.get("log.level") == "DEBUG"

    # The following are default values from nplinker_default.toml
    assert config.get("log.file") is None
    assert config.log.use_console is True

    assert config.mibig.to_use is True
    assert config.mibig.version == "3.1"

    assert (
        config.bigscape.parameters
        == "--mibig --clans-off --mix --include_singletons --cutoffs 0.30"
    )
    assert config.bigscape.cutoff == "0.30"
    assert config.bigscape.version == 1

    assert config.scoring.methods == ["metcalf"]
