import json
import os
from pathlib import Path
import pytest
from nplinker.strain import Strain
from nplinker.strain import StrainCollection
from nplinker.strain.utils import create_strain_mappings
from nplinker.strain.utils import extract_features_metabolome_id
from nplinker.strain.utils import extract_strain_metadata
from nplinker.strain.utils import load_user_strains
from nplinker.strain.utils import merge_bgcs_features
from nplinker.strain.utils import podp_generate_strain_mappings


@pytest.fixture
def user_strains_file(tmp_path):
    """Create a JSON file containing user specified strains."""
    data = {
        "strain_ids": ["strain1", "strain2", "strain3"],
    }
    file_path = tmp_path / "user_strains.json"
    with open(file_path, "w") as f:
        json.dump(data, f)
    return file_path


def test_load_user_strains(user_strains_file):
    """Test load_user_strains function."""
    actual = load_user_strains(user_strains_file)
    expected = {Strain("strain1"), Strain("strain2"), Strain("strain3")}
    assert actual == expected


def test_podp_generate_strain_mappings(monkeypatch, tmp_path):
    # mock functions called by the tested function
    mappings_strain_bgc = {
        "strain1": {"bgc1", "bgc2"},
        "strain2": {"bgc3"},
    }
    mappings_strain_spectrum = {"strain1": {"spec1", "spec2"}, "strain2": {"spec3"}}

    # monkeypatch requires the mocked function is in the same scope of the tested function
    monkeypatch.setattr(
        "nplinker.strain.utils.extract_mappings_strain_id_original_genome_id",
        lambda *args: {},
    )  # any return value is fine
    monkeypatch.setattr(
        "nplinker.strain.utils.extract_mappings_original_genome_id_resolved_genome_id",
        lambda *args: {},
    )
    monkeypatch.setattr(
        "nplinker.strain.utils.extract_mappings_resolved_genome_id_bgc_id",
        lambda *args: {},
    )
    monkeypatch.setattr(
        "nplinker.strain.utils.get_mappings_strain_id_bgc_id",
        lambda *args: mappings_strain_bgc,
    )

    monkeypatch.setattr(
        "nplinker.strain.utils.extract_mappings_strain_id_ms_filename",
        lambda *args: {},
    )
    monkeypatch.setattr(
        "nplinker.strain.utils.extract_mappings_ms_filename_spectrum_id",
        lambda *args: {},
    )
    monkeypatch.setattr(
        "nplinker.strain.utils.get_mappings_strain_id_spectrum_id",
        lambda *args: mappings_strain_spectrum,
    )

    # Create the expected
    expected_dict = {"strain1": {"bgc1", "bgc2", "spec1", "spec2"}, "strain2": {"bgc3", "spec3"}}
    expected_sc = StrainCollection()
    for strain_id, ids in expected_dict.items():
        strain = Strain(strain_id)
        for iid in ids:
            strain.add_alias(iid)
        expected_sc.add(strain)

    # Call function to generate strain mappings
    output_file = tmp_path / "output.json"
    result = podp_generate_strain_mappings(
        "dummy_podp_project_file",
        "dummy_genome_status_file",
        "dummy_genome_bgc_mappings_file",
        "dummy_gnps_file_mapping_file",
        output_file,
    )
    # check returned value
    assert isinstance(result, StrainCollection)
    assert result == expected_sc
    # check output file
    sc = StrainCollection.read_json(output_file)
    assert sc == expected_sc

class FakeAntismashBGCLoader:
    def __init__(self, bgc_path):
        pass
    def get_genome_bgcs_mapping(self):
        # create the bgc_dict for testing
        return {
            "genome_1": ["bgc_1", "bgc_2"],
            "genome_2": ["bgc_3"],
            "genome_3": ["bgc_4", "bgc_5"],
        }

# redefine function to use  FakeAntismashBGCLoader instead
def extract_bgcs_genome_id_test(strain_genome, bgc_path):
    """Extract bgcs based on the strain_genome mapping.

    Args:
        strain_genome: dict that comes from extract_strain_metadata function
        bgc_path: path of the folder of antismash results
    """
    bgc_loader = FakeAntismashBGCLoader(bgc_path)
    bgc_dict = bgc_loader.get_genome_bgcs_mapping()

    strain_bgcs = {}

    for strain_id, genome_id in strain_genome.items():
        if genome_id in bgc_dict:
            strain_bgcs[strain_id] = bgc_dict[genome_id]

    return strain_bgcs

def test_extract_bgcs_genome_id():
    strain_genome = {
        "strain_1": "genome_1",
        "strain_2": "genome_2",
        "strain_3": "genome_3",
        "strain_4": "genome_4",
    }
    # now test everything together
    strain_bgcs = extract_bgcs_genome_id_test(strain_genome, None)  # None to avoid path access

    expected_strain_bgcs = {
        "strain_1": ["bgc_1", "bgc_2"],
        "strain_2": ["bgc_3"],
        "strain_3": ["bgc_4", "bgc_5"],
    }
    assert strain_bgcs == expected_strain_bgcs

def test_extract_strain_metadata(tmp_path):
    # creation of a tab-separated file for genome
    data_genome = """StrainID\tGenomeID
    strain1\tgenome1
    strain2\tgenome2"""
    file_path = tmp_path / "strain_metadata_genomics.txt"
    assert file_path.suffix == ".txt", "Input file should have a .txt extension"
    with open(file_path, "w") as f:
        f.write(data_genome)

    actual_genome = extract_strain_metadata(file_path)
    assert isinstance(actual_genome, dict), "The output should be a dictionary"
    assert all(
        isinstance(k, str) for k in actual_genome.keys()
    ), "All keys in the output dictionary should be StrainIDs (strings)"
    expected_genome = {
        "strain1": ["genome1"],
        "strain2": ["genome2"],
    }
    assert actual_genome == expected_genome

    # creation of a tab-separated file for spectra
    data_spectra = """StrainID\tSpectraID
    strain1\tspectra1.mzML
    strain2\tspectra2.mzML"""
    file_path = tmp_path / "strain_metadata_metabolomics.txt"
    assert file_path.suffix == ".txt", "Input file should have a .txt extension"

    with open(file_path, "w") as f:
        f.write(data_spectra)

    actual_spectra = extract_strain_metadata(file_path)
    assert isinstance(actual_spectra, dict), "The output should be a dictionary"
    assert all(
        isinstance(k, str) for k in actual_spectra.keys()
    ), "All keys in the output dictionary should be StrainIDs (strings)"
    expected_spectra = {
        "strain1": ["spectra1.mzML"],
        "strain2": ["spectra2.mzML"],
    }
    assert actual_spectra == expected_spectra

def extract_mappings_ms_filename_spectrum_id(_):
    # simulating the output function
    return {
        "spectrum1": ["featureA", "featureB"],
        "spectrum2": ["featureC"],
        "spectrum3": ["featureA", "featureD"],
    }


def extract_features_metabolome_id_test(strain_spectra, _):
    """Fake function"""
    features_dict = extract_mappings_ms_filename_spectrum_id(
        None
    )  # None or ignore it completely!!!
    strain_features = {}
    for strain_id, spectra in strain_spectra.items():
        if strain_id == "StrainID":
            continue
        if isinstance(spectra, str):
            spectra = [spectra]
        features_set = set()

        for spectrum in spectra:
            if spectrum in features_dict:
                features_set.update(features_dict[spectrum])

        strain_features[strain_id] = sorted(features_set)
    return strain_features


def test_extract_features_metabolome_id():
    # Step 1: Prepare the strain_spectra data
    strain_spectra = {
        "StrainID": "ExtractID",  # This is ignored
        "Strain1": ["15b.mzXML", "12c.mzXML"],
        "Strain2": "15a.mzXML"
    }

    # Get the absolute path to the test file
    test_file = Path(__file__).parent.parent / "data/gnps/nodes.tsv"

    # Call the function with the dynamically determined path
    result = extract_features_metabolome_id(strain_spectra, str(test_file))

    # Check if the result matches the expected output
    assert len(result) == 2



def test_merge_bgcs_features():
    strain_bgcs_fake= {'strain_1': ['bgc1','bgc2','bgc3','bgc4',]}
    strain_features_fake= {'strain_1': ['feature1','feature2','feature3','feature4']}
    strain_bgcs_features = merge_bgcs_features(strain_bgcs_fake, strain_features_fake)
    expected = {'strain_1': ['bgc1','bgc2','bgc3','bgc4','feature1','feature2','feature3','feature4']}
    assert strain_bgcs_features == expected, f"Test failed! Expected {expected}, but got {strain_bgcs_features}"

def test_create_strain_mappings():
    strain_bgcs = {
        "Strain1": ["BGC1", "BGC2"],
        "Strain2": ["BGC3"],
    }

    strain_features = {"Strain1": ["featureA", "featureB"], "Strain2": ["featureC"]}

    expected_output = {
        "version": "1.0",
        "strain_mappings": [
            {"strain_id": "Strain1", "strain_alias": ["BGC1", "BGC2", "featureA", "featureB"]},
            {"strain_id": "Strain2", "strain_alias": ["BGC3", "featureC"]},
        ],
    }
    result = create_strain_mappings(
        strain_bgcs, strain_features, version="1.0", filename="strain_mappings.json"
    )
    assert result == expected_output

    # check if the json file was created
    file_path = "strain_mappings.json"
    assert os.path.exists(file_path), f"File {file_path} was not created."
    with open(file_path, "r") as json_file:
        file_content = json.load(json_file)
        assert (
            file_content == expected_output
        ), f"File content {file_content} does not match expected output {expected_output}"
