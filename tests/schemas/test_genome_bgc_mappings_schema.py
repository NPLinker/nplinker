from jsonschema import validate
from jsonschema.exceptions import ValidationError
import pytest
from nplinker.schemas import GENOME_BGC_MAPPINGS_SCHEMA


# Note:
# The function `validate` will first verify that the provided schema is itself valid (SchemaError exception),
# and then check that the given data is valid against the schema (ValidationError exception).
# It's assumed that the schema is valid, so we only need to test the ValidationError exception.
# see https://python-jsonschema.readthedocs.io/en/stable/validate/#the-basics

# Prepare invalid data
data_no_mappings = {"version": "1.0"}

data_empty_mappings = {"mappings": [], "version": "1.0"}

data_no_genome_id = {"mappings": [{"BGC_ID": ["bgc1", "bgc2"]}], "version": "1.0"}

data_empty_genome_id = {"mappings": [{"genome_ID": "", "BGC_ID": ["bgc1"]}], "version": "1.0"}

data_invalid_genome_id = {"mappings": [{"genome_ID": 1, "BGC_ID": ["bgc1"]}], "version": "1.0"}

data_no_bgc_id = {"mappings": [{"genome_ID": "genome1"}], "version": "1.0"}

data_empty_bgc_id_list = {"mappings": [{"genome_ID": "genome1", "BGC_ID": []}], "version": "1.0"}

data_empty_bgc_id = {"mappings": [{"genome_ID": "genome1", "BGC_ID": [""]}], "version": "1.0"}

data_invalid_bgc_id = {"mappings": [{"genome_ID": "genome1", "BGC_ID": [1]}], "version": "1.0"}

data_duplicate_bgc_id = {
    "mappings": [{"genome_ID": "genome1", "BGC_ID": ["bgc1", "bgc1"]}],
    "version": "1.0",
}

data_no_version = {"mappings": [{"genome_ID": "genome1", "BGC_ID": ["bgc1", "bgc2"]}]}

data_empty_version = {
    "mappings": [{"genome_ID": "genome1", "BGC_ID": ["bgc1", "bgc2"]}],
    "version": "",
}

data_invalid_version = {
    "mappings": [{"genome_ID": "genome1", "BGC_ID": ["bgc1"]}],
    "version": "1.0.0",
}


# Test schema aginast invalid data
@pytest.mark.parametrize(
    "data, expected",
    [
        [data_no_mappings, "'mappings' is a required property"],
        [data_empty_mappings, "[] is too short"],
        [data_no_genome_id, "'genome_ID' is a required property"],
        [data_empty_genome_id, "'' is too short"],
        [data_invalid_genome_id, "1 is not of type 'string'"],
        [data_no_bgc_id, "'BGC_ID' is a required property"],
        [data_empty_bgc_id_list, "[] is too short"],
        [data_empty_bgc_id, "'' is too short"],
        [data_invalid_bgc_id, "1 is not of type 'string'"],
        [data_duplicate_bgc_id, "['bgc1', 'bgc1'] has non-unique elements"],
        [data_no_version, "'version' is a required property"],
        [data_empty_version, "'' is not one of ['1.0']"],
        [data_invalid_version, "'1.0.0' is not one of ['1.0']"],
    ],
)
def test_invalid_data(data, expected):
    with pytest.raises(ValidationError) as e:
        validate(data, GENOME_BGC_MAPPINGS_SCHEMA)
    assert e.value.message == expected


# Test schema aginast valid data
def test_valid_data():
    data = {
        "mappings": [
            {"genome_ID": "genome1", "BGC_ID": ["bgc1", "bgc2"]},
            {"genome_ID": "genome2", "BGC_ID": ["bgc3"]},
        ],
        "version": "1.0",
    }
    try:
        validate(data, GENOME_BGC_MAPPINGS_SCHEMA)
    except ValidationError:
        pytest.fail("Unexpected ValidationError")
