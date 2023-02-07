from pathlib import Path
import pytest
from nplinker.loader import DatasetLoader
from . import DATA_DIR

@pytest.fixture
def config():
    return {
        "dataset" : {
            "root": DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra",
            "platform_id": ""
        }
    }

def test_default(config):
    sut = DatasetLoader(config)
    assert sut._platform_id == config["dataset"]["platform_id"]


def test_has_spectra_path(config):
    sut = DatasetLoader(config)
    sut._init_metabolomics_paths()
    actual = sut.mgf_file
    assert actual == str(config["dataset"]["root"] / "METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf")
