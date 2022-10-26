import os
from nplinker.metabolomics.GNPSMolecularFamilyLoader import GNPSMolecularFamilyLoader
from . import DATA_DIR


def test_has_molecular_families():
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    sut = GNPSMolecularFamilyLoader(filename)
    actual = sut.families()
    assert len(actual[0]) == 19
