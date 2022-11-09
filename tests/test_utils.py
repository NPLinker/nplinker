import pytest
from nplinker.utils import find_delimiter

from tests import DATA_DIR

@pytest.mark.parametrize("filename, expected", [
    [DATA_DIR / "nodes.tsv", '\t'],
    [DATA_DIR / "nodes_mwe.csv", ',']
])
def test_find_delimiter(filename, expected):
    actual = find_delimiter(filename)
    assert actual == expected

