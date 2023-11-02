import pytest
from jsonschema import validate
from jsonschema.exceptions import ValidationError
from nplinker.schemas import GENOME_STATUS_SCHEMA


# Prepare invalid data
data_no_genome_status = {"version": "1.0"}

data_empty_genome_status = {"genome_status": [], "version": "1.0"}

data_no_original_id = {
    "genome_status": [
        {"resolved_refseq_id": "id1_refseq", "resolve_attempted": True, "bgc_path": ""}
    ],
    "version": "1.0",
}

data_empty_original_id = {
    "genome_status": [
        {
            "original_id": "",
            "resolved_refseq_id": "id1_refseq",
            "resolve_attempted": True,
            "bgc_path": "",
        }
    ],
    "version": "1.0",
}

data_invalid_original_id = {
    "genome_status": [
        {
            "original_id": 1,
            "resolved_refseq_id": "id1_refseq",
            "resolve_attempted": True,
            "bgc_path": "",
        }
    ],
    "version": "1.0",
}

data_no_resolved_refseq_id = {
    "genome_status": [{"original_id": "id1", "resolve_attempted": True, "bgc_path": ""}],
    "version": "1.0",
}

data_invalid_resolved_refseq_id = {
    "genome_status": [
        {"original_id": "id1", "resolved_refseq_id": 1, "resolve_attempted": True, "bgc_path": ""}
    ],
    "version": "1.0",
}

data_no_resolve_attempted = {
    "genome_status": [{"original_id": "id1", "resolved_refseq_id": "id1_refseq", "bgc_path": ""}],
    "version": "1.0",
}

data_invalid_resolve_attempted = {
    "genome_status": [
        {
            "original_id": "id1",
            "resolved_refseq_id": "id1_refseq",
            "resolve_attempted": 1,
            "bgc_path": "",
        }
    ],
    "version": "1.0",
}

data_no_bgc_path = {
    "genome_status": [
        {"original_id": "id1", "resolved_refseq_id": "id1_refseq", "resolve_attempted": True}
    ],
    "version": "1.0",
}

data_invalid_bgc_path = {
    "genome_status": [
        {
            "original_id": "id1",
            "resolved_refseq_id": "id1_refseq",
            "resolve_attempted": True,
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
    "version": "" "",
}

data_invalid_version = {
    "genome_status": [{"strain_id": "strain1", "strain_alias": ["alias1"]}],
    "version": "1.0.0",
}


# Test schema aginast invalid data
@pytest.mark.parametrize(
    "data, expected",
    [
        [data_no_genome_status, "'genome_status' is a required property"],
        [data_empty_genome_status, "[] is too short"],
        [data_no_original_id, "'original_id' is a required property"],
        [data_empty_original_id, "'' is too short"],
        [data_invalid_original_id, "1 is not of type 'string'"],
        [data_no_resolved_refseq_id, "'resolved_refseq_id' is a required property"],
        [data_invalid_resolved_refseq_id, "1 is not of type 'string'"],
        [data_no_resolve_attempted, "'resolve_attempted' is a required property"],
        [data_invalid_resolve_attempted, "1 is not of type 'boolean'"],
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


# Test schema aginast valid data
def test_valid_data():
    data = {
        "genome_status": [
            {
                "original_id": "id1",
                "resolved_refseq_id": "id1_refseq",
                "resolve_attempted": True,
                "bgc_path": "",
            },
            {
                "original_id": "id2",
                "resolved_refseq_id": "id2_refseq",
                "resolve_attempted": False,
                "bgc_path": "",
            },
        ],
        "version": "1.0",
    }
    try:
        validate(data, GENOME_STATUS_SCHEMA)
    except ValidationError:
        pytest.fail("Unexpected ValidationError")
