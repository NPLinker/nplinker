from pathlib import Path
import pytest
from nplinker.loader import DatasetLoader
from . import DATA_DIR

@pytest.fixture
def config(tmp_path: Path):
    return {
        "dataset" : {
            "root": DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra",
            "platform_id": ""
        }
    }

def test_default(config):
    sut = DatasetLoader(config)
    assert sut is not None


def test_has_spectra_path(config):
    sut = DatasetLoader(config)
    sut.validate()
    assert sut.mgf_file == config["dataset"]["root"] / "spectra.mgf"