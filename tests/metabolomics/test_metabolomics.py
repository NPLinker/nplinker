import pytest
from nplinker.metabolomics.metabolomics import make_families
from nplinker.metabolomics.metabolomics import load_dataset
from nplinker.metabolomics.metabolomics import load_spectra
from nplinker.metabolomics.molecular_family import MolecularFamily
from nplinker.metabolomics.spectrum import Spectrum
from nplinker.strain_collection import StrainCollection
from .. import DATA_DIR


@pytest.fixture
def spec_dict() -> dict[int, Spectrum]:
    mgf_file = DATA_DIR / "spectra.mgf"
    edges_file = DATA_DIR / "edges.pairsinfo"
    return load_spectra(mgf_file, edges_file)

@pytest.fixture
def spec_with_families(spec_dict) -> dict[int, Spectrum]:
    make_families(spec_dict.values())
    return spec_dict

@pytest.fixture
def molecular_families(spec_dict) -> list[MolecularFamily]:
    return make_families(spec_dict.values())


def test_load_spectra(spec_dict):
    assert len(spec_dict.keys()) > 0


def test_load_dataset():
    strains = StrainCollection()
    mgf_file = DATA_DIR / "spectra.mgf"
    edges_file = DATA_DIR /  "edges.pairsinfo"
    nodes_file = DATA_DIR / "nodes.tsv"

    spec_dict, spectra, _, _ = load_dataset(
        strains,
        mgf_file,
        edges_file,
        nodes_file
    )

    assert isinstance(spec_dict, dict)
    assert len(spectra) > 1


def test_make_families(spec_dict):
    families = make_families(spec_dict.values())
    assert len(families) == 25769
