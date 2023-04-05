import pytest
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain
from tests import DATA_DIR


@pytest.fixture
def collection(strain: Strain) -> StrainCollection:
    sut = StrainCollection()
    sut.add(strain)
    return sut


def test_repr(collection: StrainCollection):
    assert repr(collection) == str(collection)


def test_str(collection: StrainCollection):
    assert str(collection) == 'StrainCollection(n=1) [strain_1]'


def test_len(collection: StrainCollection):
    assert len(collection) == 1


def test_eq(collection: StrainCollection, strain: Strain):
    other = StrainCollection()
    other.add(strain)
    assert collection == other


def test_contains(collection: StrainCollection, strain: Strain):
    assert strain in collection
    assert strain.id in collection
    for alias in strain.aliases:
        assert alias in collection
    assert "strain_not_exist" not in collection


def test_iter(collection: StrainCollection, strain: Strain):
    for actual in collection:
        assert actual == strain

def test_add(strain: Strain):
    sut = StrainCollection()
    sut.add(strain)
    assert strain in sut
    for alias in strain.aliases:
        assert alias in sut
    assert sut._strain_dict_index[0] == strain


def test_remove(collection: StrainCollection, strain: Strain):
    assert strain in collection
    collection.remove(strain)
    with pytest.raises(KeyError):
        _ = collection._strain_dict_id[strain.id]
    assert strain not in collection
    # TODO: issue #90
    # with pytest.raises(KeyError):
    #     collection.lookup_index(0)


def test_filter(collection: StrainCollection, strain: Strain):
    assert strain in collection
    collection.add(Strain("strain_2"))
    collection.filter({strain})
    assert strain in collection
    assert "strain_2" not in collection
    assert len(collection) == 1


def test_lookup_index(collection: StrainCollection, strain: Strain):
    actual = collection.lookup_index(0)
    assert actual == strain
    with pytest.raises(KeyError):
        collection.lookup_index(1)


def test_lookup(collection: StrainCollection, strain: Strain):
    assert collection.lookup(strain.id) == strain
    for alias in strain.aliases:
        assert collection.lookup(alias) == strain
    with pytest.raises(KeyError):
        collection.lookup("strain_not_exist")


def test_add_from_file():
    sut = StrainCollection()
    sut.add_from_file(DATA_DIR / "strain_mappings.csv")
    assert len(sut) == 27
    assert len(sut.lookup_index(1).aliases) == 29


def test_save_to_file(collection: StrainCollection, tmp_path):
    collection.add(Strain("strain_2"))
    path = tmp_path / "test.csv"
    collection.save_to_file(path)
    assert path.exists()
    with open(path) as f:
        lines = f.readlines()
        assert len(lines) == 2
        assert lines[0].strip() == "strain_1,strain_1_a"
        assert lines[1].strip() == "strain_2"
