import os
import pytest
from nplinker.metabolomics.metabolomics import load_dataset
from nplinker.metabolomics.metabolomics import load_spectra
from nplinker.strain_collection import StrainCollection


testdata_dir = os.path.join(os.getcwd(),"tests", "data")


@pytest.fixture
def spec_dict() -> dict:
    mgf_file = os.path.join(testdata_dir, "spectra.mgf")
    edges_file = os.path.join(testdata_dir, "edges.pairsinfo")
    return load_spectra(mgf_file, edges_file)


def test_load_spectra(spec_dict):
    assert len(spec_dict.keys()) > 0


def test_load_dataset():
    strains = StrainCollection()
    mgf_file = os.path.join(testdata_dir, "spectra.mgf")
    edges_file = os.path.join(testdata_dir, "edges.pairsinfo")
    nodes_file = os.path.join(testdata_dir, "nodes.tsv")

    spec_dict, spectra, molfams, unknown_strains = load_dataset(
        strains,
        mgf_file,
        edges_file,
        nodes_file
    )

    assert isinstance(spec_dict, dict)
    assert len(spectra) > 1
