from nplinker.genomics import BGC
from nplinker.strains import Strain


def test_default():

    bgc = BGC(1, Strain("BGC0000001"), "BGC0000001", ["Polyketide"])
    assert bgc.id == 1
    assert isinstance(bgc.strain, Strain)
    assert bgc.strain.id == "BGC0000001"
    assert bgc.name == "BGC0000001"
    assert bgc.product_prediction == ["Polyketide"]
    assert bgc.is_mibig() == True
