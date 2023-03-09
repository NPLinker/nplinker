from pathlib import Path
import shutil
import pytest
from nplinker.metabolomics.metabolomics import load_spectra
from nplinker.metabolomics.metabolomics import make_families
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain
from nplinker.utils import extract_archive
from . import DATA_DIR


def _unpack(archive: Path):
    filepath = DATA_DIR / archive
    outdir = DATA_DIR / filepath.stem
    extract_archive(filepath, outdir)
    return filepath, outdir


@pytest.fixture(scope="session", autouse=True)
def prepare_data():
    _unpack(
        "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip"
    )
    _unpack(
        "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip"
    )
    yield
    shutil.rmtree(
        DATA_DIR /
        "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra")
    shutil.rmtree(
        DATA_DIR /
        "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data"
    )


@pytest.fixture
def spec_dict() -> dict[int, Spectrum]:
    mgf_file = DATA_DIR / "spectra.mgf"
    edges_file = DATA_DIR / "edges.pairsinfo"
    return load_spectra(mgf_file, edges_file)


@pytest.fixture
def collection_from_file() -> StrainCollection:
    filename = DATA_DIR / "strain_mappings.csv"
    sut = StrainCollection()
    sut.add_from_file(filename)
    return sut


@pytest.fixture
def strain() -> Strain:
    item = Strain("peter")
    item.add_alias("dieter")
    return item
