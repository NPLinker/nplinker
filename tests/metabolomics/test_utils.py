import pytest
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.metabolomics import add_annotation_to_spectrum
from nplinker.metabolomics import add_spectrum_to_mf
from nplinker.metabolomics import add_strains_to_spectrum
from nplinker.metabolomics import get_spectra_from_mfs
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection


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


def test_get_spectra_from_mfs(spectra):
    """Test the get_spectra_from_mfs function."""
    mf0 = MolecularFamily("mf0")
    mf0.spectra_ids = {"spec0", "spec1"}
    mf0.add_spectrum(spectra[0])
    mf0.add_spectrum(spectra[1])
    mf1 = MolecularFamily("mf1")
    mf1.spectra_ids = {
        "spec2",
        "spec-missing-1",
    }
    mf1.add_spectrum(spectra[2])
    mf2 = MolecularFamily("mf2")
    mf2.spectra_ids = {"spec-missing-2", "spec-missing-3"}
    mfs = [mf0, mf1, mf2]

    spec_from_mfs = get_spectra_from_mfs(mfs)

    assert len(spec_from_mfs) == 3
