from nplinker.config import config


def test_config_demo1():
    """Test loading the default config file (nplinker_demo1.toml)."""
    # The file "nplinker_demo1.toml" is set in ./conftest.py

    assert config.mode == "podp"
    assert config.podp_id == "4b29ddc3-26d0-40d7-80c5-44fb6631dbf9.4"

    # The following are default values from nplinker_default.toml
    assert config.log.level == "INFO"
    assert config["log.level"] == "INFO"
    assert config.get("log.level") == "INFO"
    assert config.get("log.file") is None
    assert config.log.to_stdout is True

    assert config.mibig.to_use is True
    assert config.mibig.version == "3.1"

    assert (
        config.bigscape.parameters
        == "--mibig --clans-off --mix --include_singletons --cutoffs 0.30"
    )
    assert config.bigscape.cutoff == 0.30

    assert config.scoring.methods == ["metcalf"]
