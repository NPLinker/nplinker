from jsonschema import validate
from jsonschema.exceptions import ValidationError
import pytest
from nplinker.schemas import STRAIN_MAPPINGS_SCHEMA


# Prepare invalid data
data_no_mappings = {"version": "1.0"}

data_empty_mappings = {"strain_mappings": [], "version": "1.0"}

data_no_strain_id = {"strain_mappings": [{"strain_alias": ["alias1", "alias2"]}], "version": "1.0"}

data_empty_strain_id = {
    "strain_mappings": [{"strain_id": "", "strain_alias": ["alias1"]}],
    "version": "1.0",
}

data_invalid_strain_id = {
    "strain_mappings": [{"strain_id": 1, "strain_alias": ["alias1"]}],
    "version": "1.0",
}

data_no_strain_alias = {"strain_mappings": [{"strain_id": "strain1"}], "version": "1.0"}

data_empty_strain_alias_list = {
    "strain_mappings": [{"strain_id": "strain1", "strain_alias": []}],
    "version": "1.0",
}

data_empty_strain_alias = {
    "strain_mappings": [{"strain_id": "strain1", "strain_alias": [""]}],
    "version": "1.0",
}

data_invalid_strain_alias = {
    "strain_mappings": [{"strain_id": "strain1", "strain_alias": [1]}],
    "version": "1.0",
}

data_duplicate_strain_alias = {
    "strain_mappings": [{"strain_id": "strain1", "strain_alias": ["alias1", "alias1"]}],
    "version": "1.0",
}

data_no_version = {
    "strain_mappings": [{"strain_id": "strain1", "strain_alias": ["alias1", "alias2"]}]
}

data_empty_version = {
    "strain_mappings": [{"strain_id": "strain1", "strain_alias": ["alias1", "alias2"]}],
    "version": "" "",
}

data_invalid_version = {
    "strain_mappings": [{"strain_id": "strain1", "strain_alias": ["alias1"]}],
    "version": "1.0.0",
}


# Test schema aginast invalid data
@pytest.mark.parametrize(
    "data, expected",
    [
        [data_no_mappings, "'strain_mappings' is a required property"],
        [data_empty_mappings, "[] is too short"],
        [data_no_strain_id, "'strain_id' is a required property"],
        [data_empty_strain_id, "'' is too short"],
        [data_invalid_strain_id, "1 is not of type 'string'"],
        [data_no_strain_alias, "'strain_alias' is a required property"],
        [data_empty_strain_alias_list, "[] is too short"],
        [data_empty_strain_alias, "'' is too short"],
        [data_invalid_strain_alias, "1 is not of type 'string'"],
        [data_duplicate_strain_alias, "['alias1', 'alias1'] has non-unique elements"],
        [data_no_version, "'version' is a required property"],
        [data_empty_version, "'' is not one of ['1.0']"],
        [data_invalid_version, "'1.0.0' is not one of ['1.0']"],
    ],
)
def test_invalid_data(data, expected):
    with pytest.raises(ValidationError) as e:
        validate(data, STRAIN_MAPPINGS_SCHEMA)
    assert e.value.message == expected


# Test schema aginast valid data
def test_valid_data():
    data = {
        "strain_mappings": [
            {"strain_id": "strain1", "strain_alias": ["alias1", "alias2"]},
            {"strain_id": "strain2", "strain_alias": ["alias3"]},
        ],
        "version": "1.0",
    }
    try:
        validate(data, STRAIN_MAPPINGS_SCHEMA)
    except ValidationError:
        pytest.fail("Unexpected ValidationError")
