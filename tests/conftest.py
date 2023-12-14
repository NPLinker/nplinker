import pytest
from nplinker.globals import STRAIN_MAPPINGS_FILENAME
from nplinker.metabolomics import Spectrum
from nplinker.metabolomics.metabolomics import load_spectra
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection
from . import DATA_DIR
from . import GNPS_DATA_DIR


@pytest.fixture
def spec_dict() -> dict[str, Spectrum]:
    mgf_file = GNPS_DATA_DIR / "spectra.mgf"
    edges_file = GNPS_DATA_DIR / "edges.pairsinfo"
    return load_spectra(mgf_file, edges_file)


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
