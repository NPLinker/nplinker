import pytest
from jsonschema import validate
from jsonschema.exceptions import ValidationError
from nplinker.schemas import GENOME_STATUS_SCHEMA


# Prepare invalid data
data_no_genome_status = {"version": "1.0"}

data_empty_genome_status = {"genome_status": [], "version": "1.0"}

data_no_original_id = {
    "genome_status": [{"resolved_id": "id1_refseq", "failed_previously": True, "bgc_path": ""}],
    "version": "1.0",
}

data_empty_original_id = {
    "genome_status": [
        {
            "original_id": "",
            "resolved_id": "id1_refseq",
            "failed_previously": True,
            "bgc_path": "",
        }
    ],
    "version": "1.0",
}

data_invalid_original_id = {
    "genome_status": [
        {
            "original_id": 1,
            "resolved_id": "id1_refseq",
            "failed_previously": True,
            "bgc_path": "",
        }
    ],
    "version": "1.0",
}

data_no_resolved_id = {
    "genome_status": [{"original_id": "id1", "failed_previously": True, "bgc_path": ""}],
    "version": "1.0",
}

data_invalid_resolved_id = {
    "genome_status": [
        {"original_id": "id1", "resolved_id": 1, "failed_previously": True, "bgc_path": ""}
    ],
    "version": "1.0",
}

data_no_failed_previously = {
    "genome_status": [{"original_id": "id1", "resolved_id": "id1_refseq", "bgc_path": ""}],
    "version": "1.0",
}

data_invalid_failed_previously = {
    "genome_status": [
        {
            "original_id": "id1",
            "resolved_id": "id1_refseq",
            "failed_previously": 1,
            "bgc_path": "",
        }
    ],
    "version": "1.0",
}

data_no_bgc_path = {
    "genome_status": [
        {"original_id": "id1", "resolved_id": "id1_refseq", "failed_previously": True}
    ],
    "version": "1.0",
}

data_invalid_bgc_path = {
    "genome_status": [
        {
            "original_id": "id1",
            "resolved_id": "id1_refseq",
            "failed_previously": True,
            "bgc_path": 1,
        }
    ],
    "version": "1.0",
}

data_no_version = {
    "genome_status": [{"strain_id": "strain1", "strain_alias": ["alias1", "alias2"]}]
}

data_empty_version = {
    "genome_status": [{"strain_id": "strain1", "strain_alias": ["alias1", "alias2"]}],
    "version": "",
}

data_invalid_version = {
    "genome_status": [{"strain_id": "strain1", "strain_alias": ["alias1"]}],
    "version": "1.0.0",
}


# Test schema against invalid data
@pytest.mark.parametrize(
    "data, expected",
    [
        [data_no_genome_status, "'genome_status' is a required property"],
        [data_empty_genome_status, "[] should be non-empty"],
        [data_no_original_id, "'original_id' is a required property"],
        [data_empty_original_id, "'' should be non-empty"],
        [data_invalid_original_id, "1 is not of type 'string'"],
        [data_no_resolved_id, "'resolved_id' is a required property"],
        [data_invalid_resolved_id, "1 is not of type 'string'"],
        [data_no_failed_previously, "'failed_previously' is a required property"],
        [data_invalid_failed_previously, "1 is not of type 'boolean'"],
        [data_no_bgc_path, "'bgc_path' is a required property"],
        [data_invalid_bgc_path, "1 is not of type 'string'"],
        [data_no_version, "'version' is a required property"],
        [data_empty_version, "'' is not one of ['1.0']"],
        [data_invalid_version, "'1.0.0' is not one of ['1.0']"],
    ],
)
def test_invalid_data(data, expected):
    with pytest.raises(ValidationError) as e:
        validate(data, GENOME_STATUS_SCHEMA)
    assert e.value.message == expected


# Test schema against valid data
def test_valid_data():
    data = {
        "genome_status": [
            {
                "original_id": "id1",
                "resolved_id": "id1_refseq",
                "failed_previously": True,
                "bgc_path": "",
            },
            {
                "original_id": "id2",
                "resolved_id": "id2_refseq",
                "failed_previously": False,
                "bgc_path": "",
            },
        ],
        "version": "1.0",
    }
    try:
        validate(data, GENOME_STATUS_SCHEMA)
    except ValidationError:
        pytest.fail("Unexpected ValidationError")
