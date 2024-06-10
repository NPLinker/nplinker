import pytest
from nplinker.scoring import Score


def test_init():
    s = Score("metcalf", 1.0, {})
    assert s.name == "metcalf"
    assert s.value == 1.0
    assert s.parameter == {}

    s = Score("rosetta", 1.0, {})
    assert s.name == "rosetta"
    assert s.value == 1.0
    assert s.parameter == {}

    s = Score("nplclass", 1.0, {})
    assert s.name == "nplclass"
    assert s.value == 1.0
    assert s.parameter == {}


def test_post_init():
    with pytest.raises(ValueError):
        Score("invalid", 1.0, {})


def test_getitem():
    score = Score("metcalf", 1.0, {})
    assert score["name"] == "metcalf"
    assert score["value"] == 1.0
    assert score["parameter"] == {}


def test_setitem():
    # valid values
    score = Score("metcalf", 1.0, {})
    score["name"] = "rosetta"
    score["value"] = 2.0
    score["parameter"] = {"cutoff": 0.5}
    assert score.name == "rosetta"
    assert score.value == 2.0
    assert score.parameter == {"cutoff": 0.5}

    # invalid value for name
    with pytest.raises(ValueError, match=".* is not a valid value. .*"):
        score["name"] = "invalid"
