from nplinker.strains import Strain


def test_default():
    sut = Strain("peter")
    assert sut.id == "peter"
    assert isinstance(sut.aliases, set)
    assert len(sut.aliases) == 0


def test_add_alias(strain: Strain):
    strain.add_alias("test")
    assert len(strain.aliases) == 2


def test_equal(strain: Strain):
    other = Strain("peter")
    other.add_alias("dieter")
    assert strain == other
