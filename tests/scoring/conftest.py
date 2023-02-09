import pytest
from nplinker.metabolomics.metabolomics import make_families
from nplinker.metabolomics.spectrum import Spectrum


@pytest.fixture
def spec_with_families(spec_dict) -> dict[int, Spectrum]:
    make_families(spec_dict.values())
    return spec_dict
