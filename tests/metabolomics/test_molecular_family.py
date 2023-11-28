import pytest
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection


@pytest.fixture()
def spectrum1():
    """Return a Spectrum object."""
    spec = Spectrum(spectrum_id="spec001", peaks=[(1.0, 1.0)], precursor_mz=100.0)
    spec.strains = StrainCollection()
    spec.strains.add(Strain("strain001"))
    yield spec


@pytest.fixture()
def spectrum2():
    """Return a Spectrum object."""
    spec = Spectrum(spectrum_id="spec002", peaks=[(1.0, 1.0)], precursor_mz=100.0)
    spec.strains = StrainCollection()
    spec.strains.add(Strain("strain002"))
    yield spec


def test_init():
    """Test MolecularFamily class initialization."""
    mf = MolecularFamily("mf001")
    assert mf.family_id == "mf001"
    assert mf.spectra_ids == set()
    assert mf.spectra == set()
    assert mf.strains == StrainCollection()


def test_str_repr():
    """Test __str__ and __repr__ methods."""
    mf = MolecularFamily("mf001")
    assert (
        str(mf)
        == "MolecularFamily(family_id=mf001, #Spectrum_objects=0, #spectrum_ids=0, #strains=0)"
    )
    assert repr(mf) == str(mf)


def test_eq():
    """Test __eq__ method."""
    mf1 = MolecularFamily("mf001")
    mf2 = MolecularFamily("mf001")
    assert mf1 == mf2


def test_hash():
    """Test __hash__ method."""
    mf = MolecularFamily("mf001")
    assert hash(mf) == hash("mf001")


def test_spectra(spectrum1):
    """Test spectra property."""
    mf = MolecularFamily("mf001")
    assert mf.spectra == set()
    mf.add_spectrum(spectrum1)
    assert spectrum1 in mf.spectra
    mf.spectra.remove(spectrum1)
    assert len(mf.spectra) == 0


def test_strains(spectrum1):
    """Test strains property."""
    mf = MolecularFamily("mf001")
    assert mf.strains == StrainCollection()
    mf.add_spectrum(spectrum1)
    assert Strain("strain001") in mf.strains
    mf.strains.remove(Strain("strain001"))
    assert len(mf.strains) == 0


def test_add_spectrum(spectrum1, spectrum2):
    """Test add_spectrum method."""
    mf = MolecularFamily("mf001")
    mf.add_spectrum(spectrum1)
    assert spectrum1 in mf.spectra
    assert spectrum1.spectrum_id in mf.spectra_ids
    assert Strain("strain001") in mf.strains
    mf.add_spectrum(spectrum2)
    assert spectrum2 in mf.spectra
    assert spectrum2.spectrum_id in mf.spectra_ids
    assert Strain("strain002") in mf.strains
    assert len(mf.spectra) == 2
    assert len(mf.spectra_ids) == 2
    assert len(mf.strains) == 2


def test_detach_spectrum(spectrum1, spectrum2):
    """Test detach_spectrum method."""
    mf = MolecularFamily("mf001")
    mf.add_spectrum(spectrum1)
    mf.add_spectrum(spectrum2)
    mf.detach_spectrum(spectrum1)
    assert spectrum1 not in mf.spectra
    assert spectrum1.spectrum_id not in mf.spectra_ids
    assert Strain("strain001") not in mf.strains
    assert Strain("strain002") in mf.strains
    mf.detach_spectrum(spectrum2)
    assert spectrum2 not in mf.spectra
    assert spectrum2.spectrum_id not in mf.spectra_ids
    assert Strain("strain002") not in mf.strains
    assert len(mf.spectra) == 0
    assert len(mf.spectra_ids) == 0
    assert len(mf.strains) == 0


def test_has_strain(spectrum1, spectrum2):
    """Test has_strain method."""
    mf = MolecularFamily("mf001")
    mf.add_spectrum(spectrum1)
    assert mf.has_strain(Strain("strain001"))
    mf.add_spectrum(spectrum2)
    assert mf.has_strain(Strain("strain002"))


def test_is_singleton(spectrum1, spectrum2):
    """Test is_singleton method."""
    mf = MolecularFamily("mf001")
    assert not mf.is_singleton()  # 0 spectra
    mf.add_spectrum(spectrum1)
    assert mf.is_singleton()
    mf.add_spectrum(spectrum2)
    assert not mf.is_singleton()
