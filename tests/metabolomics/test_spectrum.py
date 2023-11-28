import pytest
from nplinker.metabolomics import Spectrum


@pytest.fixture
def spectrum() -> Spectrum:
    spec = Spectrum(1, peaks=[[10, 100], [20, 150]], spectrum_id="2", precursor_mz=30, rt=100)
    return spec


def test_constructor(spectrum):
    assert spectrum is not None
