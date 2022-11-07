import pytest
from nplinker.strains import Strain


@pytest.fixture
def strain() -> Strain:
    item = Strain("peter")
    item.aliases = set(["dieter"])
    return item

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