from nplinker.genomics import BGC
from nplinker.strains import Strain


def test_default():

    bgc = BGC("BGC0000001", ["Polyketide"])
    assert bgc.bgc_id == "BGC0000001"
    assert bgc.product_prediction == ["Polyketide"]
    assert bgc.is_mibig() is True
    assert bgc.parents == set()
    assert bgc.bigscape_classes == set()
    assert bgc.strain is None

    strain = Strain("sample_strain")
    bgc.strain = strain
    assert bgc.strain == strain
