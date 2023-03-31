import pytest
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain
from tests import DATA_DIR


@pytest.fixture
def collection(strain: Strain) -> StrainCollection:
    sut = StrainCollection()
    sut.add(strain)
    return sut


def test_default():
    sut = StrainCollection()
    assert sut is not None


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


def test_add():
    sut = StrainCollection()
    item = Strain("test_id")
    item.add_alias("blub")

    sut.add(item)

    assert sut.lookup(item.id) == item
    assert sut.lookup(next(iter(item.aliases))) == item
    assert sut.lookup_index(0) == item


def test_lookup(collection: StrainCollection, strain: Strain):
    assert collection.lookup(strain.id) == strain
    for alias in strain.aliases:
        assert collection.lookup(alias) == strain
    with pytest.raises(KeyError):
        collection.lookup("strain_not_exist")


def test_lookup_index(collection: StrainCollection, strain: Strain):
    actual = collection.lookup_index(0)
    assert actual == strain


def test_lookup_index_exception(collection: StrainCollection):
    with pytest.raises(KeyError) as exc:
        collection.lookup_index(5)
    assert isinstance(exc.value, KeyError)


def test_remove(collection: StrainCollection, strain: Strain):
    assert strain in collection
    collection.remove(strain)
    with pytest.raises(KeyError):
        collection.lookup(strain.id)
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


def test_equal(collection_from_file: StrainCollection):
    other = StrainCollection()
    other.add_from_file(DATA_DIR / "strain_mappings.csv")

    assert collection_from_file == other
