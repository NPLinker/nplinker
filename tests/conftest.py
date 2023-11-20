import pytest
from nplinker.globals import STRAIN_MAPPINGS_FILENAME
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection
from . import DATA_DIR


@pytest.fixture
def collection_from_file() -> StrainCollection:
    filename = DATA_DIR / STRAIN_MAPPINGS_FILENAME
    sut = StrainCollection.read_json(filename)
    return sut


@pytest.fixture
def strain() -> Strain:
    item = Strain("strain_1")
    item.add_alias("strain_1_a")
    return item
