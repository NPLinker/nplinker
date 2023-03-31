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


def test_add_from_file(collection_from_file: StrainCollection):
    assert len(collection_from_file) == 27
    assert len(collection_from_file.lookup_index(1).aliases) == 29


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
