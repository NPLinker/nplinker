import pytest
from nplinker.strains import Strain


def test_default():
    sut = Strain("peter")
    assert sut.id == "peter"
    assert len(sut.aliases) == 0


@pytest.mark.parametrize("alias, expected", [
    ["dieter", True],
    ["ulrich", False]
])
def test_has_alias(strain: Strain, alias: str, expected: bool):
    assert strain.has_alias(alias) == expected


def test_add_alias(strain: Strain):
    strain.add_alias("test")
    assert len(strain.aliases) == 2


def test_equal(strain: Strain):
    other = Strain("peter")
    other.add_alias("dieter")

    assert strain == other