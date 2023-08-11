from nplinker.strains import Strain


def test_default():
    sut = Strain("strain_1")
    assert sut.id == "strain_1"
    assert isinstance(sut.aliases, set)
    assert len(sut.aliases) == 0


def test_repr(strain: Strain):
    assert repr(strain) == "Strain(strain_1) [1 aliases]"


def test_str(strain: Strain):
    assert str(strain) == "Strain(strain_1) [1 aliases]"


def test_eq(strain: Strain):
    other = Strain("strain_1")
    other.add_alias("strain_1_a")
    assert strain == other


def test_hash(strain: Strain):
    assert hash(strain) == hash("strain_1")


def test_names(strain: Strain):
    assert strain.names == {"strain_1", "strain_1_a"}


def test_alias(strain: Strain):
    assert len(strain.aliases) == 1
    assert "strain_1_a" in strain.aliases


def test_add_alias(strain: Strain):
    strain.add_alias("strain_1_b")
    assert len(strain.aliases) == 2
    assert "strain_1_b" in strain.aliases
