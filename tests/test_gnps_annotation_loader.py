from nplinker.metabolomics.gnps.gnps_annotation_loader import GNPSAnnotationLoader


def test_default():
    sut = GNPSAnnotationLoader()
    assert sut is not None


def test_can_open_files():
    sut = GNPSAnnotationLoader()
    
    
