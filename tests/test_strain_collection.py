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


def test_lookup(strain: Strain):
    sut = StrainCollection()
    sut.add(strain)

    assert sut.lookup(strain.id) == strain


def test_contains(collection: StrainCollection, strain: Strain):
    assert strain in collection
    assert strain.id in collection
    assert "peter" in collection
    assert "dieter" in collection
    assert "test" not in collection


def test_lookup_index(collection: StrainCollection, strain: Strain):
    actual = collection.lookup_index(0)
    assert actual == strain


def test_lookup_index_exception(collection: StrainCollection):
    with pytest.raises(KeyError) as exc:
        collection.lookup_index(5)
    assert isinstance(exc.value, KeyError)


def test_remove(collection: StrainCollection, strain: Strain):
    collection.remove(strain)

    with pytest.raises(KeyError):
        collection.lookup(strain.id)

    assert strain not in collection

    # needs fixing, see #90
    assert collection.lookup_index(0) == strain


def test_equal(collection_from_file: StrainCollection):
    other = StrainCollection()
    other.add_from_file(DATA_DIR / "strain_mappings.csv")

    assert collection_from_file == other
