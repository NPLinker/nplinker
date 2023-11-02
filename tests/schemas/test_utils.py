from jsonschema import Draft7Validator
import pytest
from nplinker.schemas import PODP_ADAPTED_SCHEMA
from nplinker.schemas import validate_podp_json


def test_podp_adapted_schema_itself():
    validator = Draft7Validator(Draft7Validator.META_SCHEMA)
    errors = validator.iter_errors(PODP_ADAPTED_SCHEMA)
    assert list(errors) == []


def test_validate_podp_json_minimum_valid_data():
    # minimum valid data, containing only required fields
    data = {
        "version": "3",
        "metabolomics": {
            "project": {"molecular_network": "01234567890123456789012345678901"},
        },
        "genomes": [
            {"genome_label": "strain1", "genome_ID": {"RefSeq_accession": "GCF_1"}},
        ],
        "genome_metabolome_links": [
            {"metabolomics_file": "ftp://example.org/001.mzXML", "genome_label": "strain1"},
        ],
    }
    try:
        validate_podp_json(data)
    except ValueError:
        pytest.fail("Unexpected ValueError")


def test_validate_podp_json_invalid_data():
    data = {}
    with pytest.raises(ValueError) as e:
        validate_podp_json(data)

    expected = """
Not match PODP adapted schema, here are the detailed error:
  - $: 'version' is a required property
  - $: 'metabolomics' is a required property
  - $: 'genomes' is a required property
  - $: 'genome_metabolome_links' is a required property
"""

    assert e.value.args[0] == expected.strip()
