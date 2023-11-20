import json
import pytest
from nplinker.strain import Strain
from nplinker.strain_loader import load_user_strains


@pytest.fixture
def user_strains_file(tmp_path):
    """Create a JSON file containing user specified strains."""
    data = {
        "strain_ids": ["strain1", "strain2", "strain3"],
    }
    file_path = tmp_path / "user_strains.json"
    with open(file_path, "w") as f:
        json.dump(data, f)
    return file_path


def test_load_user_strains(user_strains_file):
    """Test load_user_strains function."""
    actual = load_user_strains(user_strains_file)
    expected = {Strain("strain1"), Strain("strain2"), Strain("strain3")}
    assert actual == expected
