from pathlib import Path
import shutil
import pytest
from nplinker.loader import DatasetLoader
from nplinker.loader import STRAIN_MAPPINGS_FILENAME
from nplinker.metabolomics.gnps.gnps_extractor import GNPSExtractor
from nplinker.metabolomics.gnps.gnps_spectrum_loader import GNPSSpectrumLoader
from nplinker.strain_collection import StrainCollection
from . import DATA_DIR


@pytest.fixture
def config():
    return {
        "dataset" : {
            "root": DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra",
            "platform_id": "",
            "overrides": {
                "strain_mappings_file": str(DATA_DIR / STRAIN_MAPPINGS_FILENAME)
            }
        }
    }

@pytest.fixture
def config_with_new_gnps_extractor():
    GNPSExtractor(DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip", DATA_DIR / "extracted").extract()
    yield {
        "dataset" : {
            "root": DATA_DIR / "extracted",
            "platform_id": "",
            "overrides": {
                "strain_mappings_file": str(DATA_DIR / STRAIN_MAPPINGS_FILENAME)
            }
        }
    }
    shutil.rmtree(DATA_DIR / "extracted")

def test_default(config):
    sut = DatasetLoader(config)
    assert sut._platform_id == config["dataset"]["platform_id"]


def test_has_metabolomics_paths(config):
    sut = DatasetLoader(config)
    sut._init_metabolomics_paths()
    assert sut.mgf_file == str(config["dataset"]["root"] / "METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf")
    assert sut.edges_file == str(config["dataset"]["root"] / "networkedges_selfloop" / "6da5be36f5b14e878860167fa07004d6.pairsinfo")
    assert sut.nodes_file == str(config["dataset"]["root"] / "clusterinfosummarygroup_attributes_withIDs_withcomponentID" / "d69356c8e5044c2a9fef3dd2a2f991e1.tsv")
    assert sut.annotations_dir == str(config["dataset"]["root"] / "result_specnets_DB")


def test_has_metabolomics_paths_new_gnps(config_with_new_gnps_extractor):
    sut = DatasetLoader(config_with_new_gnps_extractor)
    sut._init_metabolomics_paths()
    assert sut.mgf_file == str(config_with_new_gnps_extractor["dataset"]["root"] / "spectra.mgf")
    assert sut.edges_file == str(config_with_new_gnps_extractor["dataset"]["root"] / "molecular_families.pairsinfo")
    assert sut.nodes_file == str(config_with_new_gnps_extractor["dataset"]["root"] / "file_mappings.tsv")
    assert sut.annotations_dir == str(config_with_new_gnps_extractor["dataset"]["root"])
    assert sut.annotations_config_file == str(config_with_new_gnps_extractor["dataset"]["root"] / "annotations.tsv")


def test_load_metabolomics(config):
    sut = DatasetLoader(config)
    sut._init_paths()
    sut._load_strain_mappings()
    sut._load_metabolomics()

    expected_spectra = GNPSSpectrumLoader(sut.mgf_file).spectra()

    #HH TODO: switch to different comparison as soon as strains are implemented
    assert len(sut.spectra) == len(expected_spectra)
    assert len(sut.molfams) == 429

def test_has_strain_mappings(config):
    sut = DatasetLoader(config)
    sut._init_paths()
    assert sut.strain_mappings_file == str(DATA_DIR / STRAIN_MAPPINGS_FILENAME)


def test_load_strain_mappings(config):
    sut = DatasetLoader(config)
    sut._init_paths()
    sut._load_strain_mappings()

    actual = sut.strains
    expected = StrainCollection.read_json(sut.strain_mappings_file)

    assert actual == expected
