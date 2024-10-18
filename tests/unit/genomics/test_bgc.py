from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.strain import Strain


def test_default():
    bgc = BGC("BGC0000001", "Polyketide")
    assert bgc.id == "BGC0000001"
    assert bgc.product_prediction == ("Polyketide",)
    assert bgc.is_mibig() is True
    assert bgc.parents == set()
    assert bgc.bigscape_classes == set()
    assert bgc.strain is None

    strain = Strain("sample_strain")
    bgc.strain = strain
    assert bgc.strain == strain


def test_add_and_detach_parent():
    bgc = BGC("BGC0000001", "Polyketide")
    gcf = GCF("1")
    bgc.add_parent(gcf)
    assert bgc.parents == {gcf}
    bgc.detach_parent(gcf)
    assert bgc.parents == set()


def test_to_dict():
    bgc = BGC("BGC0000001", "Polyketide", "NRP")
    bgc.strain = Strain("sample_strain")
    bgc.description = "Sample description"

    dict_repr = bgc.to_dict()
    assert dict_repr["GCF_id"] == set()
    assert dict_repr["GCF_bigscape_class"] == set()
    assert dict_repr["BGC_name"] == "BGC0000001"
    assert dict_repr["product_prediction"] == ("Polyketide", "NRP")
    assert dict_repr["mibig_bgc_class"] is None
    assert dict_repr["description"] == "Sample description"
    assert dict_repr["strain_id"] == "sample_strain"
    assert dict_repr["antismash_id"] is None
    assert dict_repr["antismash_region"] is None

    bgc.add_parent(GCF("1"))
    bgc.mibig_bgc_class = ("NRP",)
    bgc.antismash_id = "ABC_0001"
    bgc.antismash_region = 1
    dict_repr = bgc.to_dict()
    assert dict_repr["GCF_id"] == {"1"}
    assert dict_repr["GCF_bigscape_class"] == set()
    assert dict_repr["mibig_bgc_class"] == ("NRP",)
    assert dict_repr["antismash_id"] == "ABC_0001"
    assert dict_repr["antismash_region"] == 1
