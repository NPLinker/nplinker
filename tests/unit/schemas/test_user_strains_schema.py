import pytest
from jsonschema import validate
from jsonschema.exceptions import ValidationError
from nplinker.schemas import USER_STRAINS_SCHEMA


# Test schema against invalid data
data_no_strain_ids = {"version": "1.0"}
data_empty_strain_ids = {"strain_ids": [], "version": "1.0"}
data_invalid_strain_ids = {
    "strain_ids": [
        1,
    ],
    "version": "1.0",
}
data_empty_version = {"strain_ids": ["strain1", "strain2"], "version": ""}
data_invalid_version = {"strain_ids": ["strain1", "strain2"], "version": "1.0.0"}


@pytest.mark.parametrize(
    "data, expected",
    [
        [data_no_strain_ids, "'strain_ids' is a required property"],
        [data_empty_strain_ids, "[] should be non-empty"],
        [data_invalid_strain_ids, "1 is not of type 'string'"],
        [data_empty_version, "'' is not one of ['1.0']"],
        [data_invalid_version, "'1.0.0' is not one of ['1.0']"],
    ],
)
def test_invalid_data(data, expected):
    """Test user strains schema against invalid data."""
    with pytest.raises(ValidationError) as e:
        validate(data, USER_STRAINS_SCHEMA)
    assert e.value.message == expected


# Test schema against valid data
data = {"strain_ids": ["strain1", "strain2"], "version": "1.0"}
data_no_version = {"strain_ids": ["strain1", "strain2"]}


@pytest.mark.parametrize(
    "data",
    [
        data,
        data_no_version,
    ],
)
def test_valid_data(data):
    """Test user strains schema against valid data."""
    try:
        validate(data, USER_STRAINS_SCHEMA)
    except ValidationError:
        pytest.fail("Unexpected ValidationError")
