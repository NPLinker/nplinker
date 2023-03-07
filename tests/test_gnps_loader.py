from nplinker.metabolomics.gnps.gnps_loader import GNPSLoader


def test_default():
    sut = GNPSLoader()
    assert sut is not None