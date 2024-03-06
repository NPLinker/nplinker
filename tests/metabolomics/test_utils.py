import json
import pytest
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.metabolomics.utils import add_annotation_to_spectrum
from nplinker.metabolomics.utils import add_spectrum_to_mf
from nplinker.metabolomics.utils import add_strains_to_spectrum
from nplinker.metabolomics.utils import extract_mappings_ms_filename_spectrum_id
from nplinker.metabolomics.utils import extract_mappings_strain_id_ms_filename
from nplinker.metabolomics.utils import get_mappings_strain_id_spectrum_id
from nplinker.strain import Strain
from nplinker.strain import StrainCollection


@pytest.fixture
def spectra():
    """Fixture for a list of Spectrum objects."""
    # The order of the spectra is important for the tests.
    return [
        Spectrum("spec0", [100, 200], [0.1, 0.2], 150),
        Spectrum("spec1", [100, 200], [0.1, 0.2], 150),
        Spectrum("spec2", [100, 200], [0.1, 0.2], 150),
    ]


def test_add_annotation_to_spectrum(spectra):
    """Test the add_annotation_to_spectrum function."""
    annotations = {
        "spec0": {"annotation": "annotation_0"},
        "spec1": {"annotation": "annotation_1"},
        "spec3": {"annotation": "annotation_3"},
    }

    add_annotation_to_spectrum(annotations, spectra)

    for i, spec in enumerate(spectra):
        if i < 2:
            assert spec.gnps_annotations == {"annotation": f"annotation_{i}"}
        else:
            assert spec.gnps_annotations == {}


def test_add_strains_to_spectrum(spectra):
    """Test the add_strains_to_spectrum function."""
    strains = StrainCollection()
    strain0 = Strain("spec0")  # spectrum id as strain id
    strain1 = Strain("strain1")
    strain1.add_alias("spec1")  # spectrum id as strain alias
    strains.add(strain0)
    strains.add(strain1)

    spectra_with_strains, spectra_without_strains = add_strains_to_spectrum(strains, spectra)

    assert len(spectra_with_strains) == 2
    assert len(spectra_without_strains) == 1
    assert spectra_with_strains == [spectra[0], spectra[1]]
    assert spectra_without_strains == [spectra[2]]
    assert strain0 in spectra_with_strains[0].strains
    assert strain1 in spectra_with_strains[1].strains
    assert spectra_without_strains[0].strains == StrainCollection()


def test_add_spectrum_to_mf(spectra):
    """Test the add_spectrum_to_mf function."""
    # Prepare the molecular families
    mf0 = MolecularFamily("mf0")
    mf0.spectra_ids = {"spec0", "spec1"}
    mf1 = MolecularFamily("mf1")
    mf1.spectra_ids = {
        "spec2",
        "spec-missing-1",
    }
    mf2 = MolecularFamily("mf2")
    mf2.spectra_ids = {"spec-missing-2", "spec-missing-3"}
    mfs = [mf0, mf1, mf2]

    mf_with_spec, mf_without_spec, mf_missing_spec = add_spectrum_to_mf(spectra, mfs)

    assert len(mf_with_spec) == 2
    assert len(mf_without_spec) == 1
    assert len(mf_missing_spec) == 2
    assert mf_with_spec == [mf0, mf1]
    assert mf_without_spec == [mf2]
    assert mf_missing_spec == {mf1: {"spec-missing-1"}, mf2: {"spec-missing-2", "spec-missing-3"}}
    assert mf0.spectra == {spectra[0], spectra[1]}
    assert mf1.spectra == {spectra[2]}
    assert mf2.spectra == set()


def test_extract_mappings_strain_id_ms_filename(tmp_path):
    test_data = {
        "genome_metabolome_links": [
            {"genome_label": "strain1", "metabolomics_file": "http://example.com/file1.mzXML"},
            {"genome_label": "strain1", "metabolomics_file": "http://example.com/file2.mzXML"},
            {"genome_label": "strain2", "metabolomics_file": "http://example.com/file3.mzXML"},
            {"genome_label": "strain3", "metabolomics_file": "http://example.com/file4.mzXML"},
        ],
        "genomes": [
            {"genome_label": "strain1", "genome_ID": {"RefSeq_accession": "id1"}},
        ],
        "metabolomics": {"project": {"molecular_network": "01234567890123456789012345678901"}},
        "version": "3",
    }
    test_file = tmp_path / "test_data.json"
    with open(test_file, "w") as f:
        json.dump(test_data, f)
    expected_result = {
        "strain1": {"file1.mzXML", "file2.mzXML"},
        "strain2": {"file3.mzXML"},
        "strain3": {"file4.mzXML"},
    }

    assert extract_mappings_strain_id_ms_filename(test_file) == expected_result


def test_extract_mappings_ms_filename_spectrum_id(tmp_path):
    test_data = "cluster index\tAllFiles\nspec1\tfile1.mzXML:123###\nspec2\tfile2.mzXML:123###\nspec3\tfile2.mzXML:123###file3.mzXML:123###\n"
    test_file = tmp_path / "test_data.tsv"
    with open(test_file, "w") as f:
        f.write(test_data)
    expected_result = {
        "file1.mzXML": {"spec1"},
        "file2.mzXML": {"spec2", "spec3"},
        "file3.mzXML": {"spec3"},
    }

    assert extract_mappings_ms_filename_spectrum_id(test_file) == expected_result


def test_get_mappings_strain_id_spectrum_id():
    mappings_strain_id_ms_filename = {
        "strain1": {"file1.mzXML", "file2.mzXML"},
        "strain2": {"file3.mzXML"},
        "strain3": {"file4.mzXML"},
    }
    mappings_ms_filename_spectrum_id = {
        "file1.mzXML": {"spec1"},
        "file2.mzXML": {"spec2", "spec3"},
        "file3.mzXML": {"spec3"},
    }

    expected_mappings_dict = {
        "strain1": {"spec1", "spec2", "spec3"},
        "strain2": {"spec3"},
    }
    actual_mappings_dict = get_mappings_strain_id_spectrum_id(
        mappings_strain_id_ms_filename, mappings_ms_filename_spectrum_id
    )

    assert actual_mappings_dict == expected_mappings_dict
