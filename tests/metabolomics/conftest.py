import pytest
from nplinker.metabolomics.metabolomics import make_families
from nplinker.metabolomics.molecular_family import MolecularFamily


@pytest.fixture
def molecular_families(spec_dict) -> list[MolecularFamily]:
    return make_families(spec_dict.values())
