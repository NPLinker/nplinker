from nplinker.genomics import BGC
from nplinker.genomics import GCF
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


def test_add_and_detach_parent():
    bgc = BGC("BGC0000001", ["Polyketide"])
    gcf = GCF("1")
    bgc.add_parent(gcf)
    assert bgc.parents == {gcf}
    bgc.detach_parent(gcf)
    assert bgc.parents == set()
