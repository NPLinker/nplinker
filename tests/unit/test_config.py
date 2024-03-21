from nplinker.config import config


def test_config():
    """Test loading the default config file."""
    # The default config file is set in "./conftest.py", which is "data/nplinker_local_mode.toml"
    assert config.mode == "local"
    assert config.log.level == "DEBUG"
    assert config["log.level"] == "DEBUG"
    assert config.get("log.level") == "DEBUG"

    # The following are default values from nplinker_default.toml
    assert config.get("log.file") is None
    assert config.log.to_stdout is True

    assert config.mibig.to_use is True
    assert config.mibig.version == "3.1"

    assert (
        config.bigscape.parameters
        == "--mibig --clans-off --mix --include_singletons --cutoffs 0.30"
    )
    assert config.bigscape.cutoff == "0.30"

    assert config.scoring.methods == ["metcalf"]
