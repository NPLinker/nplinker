import os

import pytest
from nplinker.metabolomics.gnps.gnps_molecular_family_loader import GNPSMolecularFamilyLoader
from . import DATA_DIR

@pytest.fixture
def molecular_families_gnps():
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    sut = GNPSMolecularFamilyLoader(filename)
    return sut.families()


def test_has_molecular_families():
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    sut = GNPSMolecularFamilyLoader(filename)
    actual = sut.families()
    assert len(actual) == 25769
    assert len(actual[0].spectra_ids) == 19
