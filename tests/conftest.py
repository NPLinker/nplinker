from pathlib import Path
import shutil
import pytest

from nplinker.utils import extract_archive

from . import DATA_DIR

def _unpack(archive: Path):
    filepath = DATA_DIR / archive
    outdir = DATA_DIR / filepath.stem
    extract_archive(filepath, outdir)
    return filepath, outdir

@pytest.fixture(scope="session", autouse=True)
def prepare_data():
    _unpack("ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip")
    _unpack("ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip")
    yield
    shutil.rmtree(DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra")
    shutil.rmtree(DATA_DIR / "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data")