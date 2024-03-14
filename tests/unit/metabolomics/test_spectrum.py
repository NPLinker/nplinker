import numpy as np
import pytest
from nplinker.metabolomics import Spectrum
from nplinker.strain import Strain
from nplinker.strain import StrainCollection


@pytest.mark.parametrize(
    "rt, metadata, expected_metadata",
    [
        [0, None, {}],
        [1, {"info": "test"}, {"info": "test"}],
    ],
)
def test_init(rt, metadata, expected_metadata):
    """Test the initialization of the Spectrum class."""
    spec = Spectrum("spec1", [100, 200], [0.1, 0.2], 150, rt, metadata)

    assert spec.spectrum_id == "spec1"
    assert spec.mz == [100, 200]
    assert spec.intensity == [0.1, 0.2]
    assert spec.precursor_mz == 150
    assert spec.rt == rt
    assert spec.metadata == expected_metadata

    # test the default values of the attributes
    assert spec.gnps_annotations == {}
    assert spec.gnps_id is None
    assert spec.strains == StrainCollection()
    assert spec.family is None


def test_str_repr():
    """Test the __str__ and __repr__ methods."""
    spec = Spectrum("spec1", [100, 200], [0.1, 0.2], 150)
    assert str(spec) == "Spectrum(spectrum_id=spec1, #strains=0)"
    assert repr(spec) == "Spectrum(spectrum_id=spec1, #strains=0)"


def test_eq():
    """Test the __eq__ method."""
    spec1 = Spectrum("spec1", [100, 200], [0.1, 0.2], 150, 0, {"info": "test"})
    spec2 = Spectrum("spec1", [100, 200], [0.1, 0.2], 150, 0, {"info": "test"})
    spec3 = Spectrum("spec2", [100, 200], [0.1, 0.2], 150, 0, {"info": "test"})

    assert spec1 == spec2
    assert spec1 != spec3


def test_hash():
    """Test the __hash__ method."""
    spec = Spectrum("spec1", [100, 200], [0.1, 0.2], 150)
    assert hash(spec) == hash(("spec1", 150))


def test_peaks():
    """Test the peaks attribute."""
    spec = Spectrum("spec1", [100, 200], [0.1, 0.2], 150)
    assert np.array_equal(spec.peaks, np.array([[100, 0.1], [200, 0.2]]))


def test_has_strain():
    """Test the has_strain method."""
    spec = Spectrum("spec1", [100, 200], [0.1, 0.2], 150)
    strain1 = Strain("strain1")
    strain2 = Strain("strain2")

    spec.strains.add(strain1)
    assert spec.has_strain(strain1)
    assert not spec.has_strain(strain2)
