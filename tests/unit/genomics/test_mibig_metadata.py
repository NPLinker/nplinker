import pytest
from nplinker.genomics.mibig import MibigMetadata
from .. import DATA_DIR


@pytest.mark.parametrize(
    "version, expected", [["v1.4", "Polyketide"], ["v3.1", "Polyketide"], ["v4.0", "PKS"]]
)
def test_versions(version, expected):
    json_file = DATA_DIR / "mibig" / f"BGC0000001_{version}.json"
    metadata = MibigMetadata(str(json_file))
    assert metadata.file == str(json_file)
    assert isinstance(metadata.metadata, dict)
    assert metadata.mibig_accession == "BGC0000001"
    assert metadata.biosyn_class == (expected,)
