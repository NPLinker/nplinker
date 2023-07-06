import json
import pytest
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain


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


def test_remove(collection: StrainCollection, strain: Strain):
    assert strain in collection
    collection.remove(strain)
    with pytest.raises(KeyError):
        _ = collection._strain_dict_name[strain.id]
    assert strain not in collection


def test_filter(collection: StrainCollection, strain: Strain):
    assert strain in collection
    collection.add(Strain("strain_2"))
    collection.filter({strain})
    assert strain in collection
    assert "strain_2" not in collection
    assert len(collection) == 1


def test_lookup(collection: StrainCollection, strain: Strain):
    assert collection.lookup(strain.id) == strain
    for alias in strain.aliases:
        assert collection.lookup(alias) == strain
    with pytest.raises(KeyError):
        collection.lookup("strain_not_exist")


@pytest.fixture
def json_file(tmp_path):
    data = {
        "strain_mappings": [{
            "strain_id": "strain_1",
            "strain_alias": ["alias_1", "alias_2"]
        }, {
            "strain_id": "strain_2",
            "strain_alias": ["alias_3", "alias_4"]
        }]
    }
    file_path = tmp_path / "test.json"
    with open(file_path, "w") as f:
        json.dump(data, f)
    return file_path


def test_read_json(json_file):
    expected_strain_1 = Strain("strain_1")
    expected_strain_1.add_alias("alias_1")
    expected_strain_1.add_alias("alias_2")
    expected_strain_2 = Strain("strain_2")
    expected_strain_2.add_alias("alias_3")
    expected_strain_2.add_alias("alias_4")
    expected_collection = StrainCollection()
    expected_collection.add(expected_strain_1)
    expected_collection.add(expected_strain_2)

    actual_collection = StrainCollection.read_json(json_file)
    assert actual_collection == expected_collection


def test_to_json(collection: StrainCollection, tmp_path):
    # tests writing to string
    expected_data = {
        "strain_mappings": [{
            "strain_id": "strain_1",
            "strain_alias": ["strain_1_a"]
        }],
        "version":
        "1.0"
    }
    expected_json = json.dumps(expected_data)
    actual_json = collection.to_json()
    assert actual_json == expected_json

    # tests writing to file
    file_path = tmp_path / "test.json"
    collection.to_json(file_path)
    with open(file_path, "r") as f:
        actual_data = json.load(f)
    assert actual_data == expected_data
