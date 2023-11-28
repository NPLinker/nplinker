import pytest
from nplinker.metabolomics import Spectrum


@pytest.fixture
def spectrum() -> Spectrum:
    spec = Spectrum(spectrum_id="2", mz=[10, 20], intensity=[100, 150], precursor_mz=30, rt=100)
    return spec


def test_constructor(spectrum):
    assert spectrum is not None
