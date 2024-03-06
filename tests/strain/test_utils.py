import json
import pytest
from nplinker.strain import Strain
from nplinker.strain import StrainCollection
from nplinker.strain.utils import load_user_strains
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
