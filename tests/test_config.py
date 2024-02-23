from dynaconf import Dynaconf
from nplinker.config import config
from . import DATA_DIR


def test_config_demo1():
    """Test loading the default config file (nplinker_demo1.toml)."""
    # the file "nplinker_demo1.toml" is set in ./conftest.py
    assert config.loglevel == "DEBUG"
    assert config["loglevel"] == "DEBUG"
    assert config.get("loglevel") == "DEBUG"

    assert config.logfile == ""
    assert config.repro_file == ""

    assert config.dataset.root == "platform:MSV000079284"
    assert config["dataset.root"] == "platform:MSV000079284"
    assert config.get("dataset.root") == "platform:MSV000079284"

    assert config.dataset.platform_id == "MSV000079284"
    assert config.webapp.tables_metcalf_score == 3.0


