from nplinker.genomics.mibig import MibigBGC


def test_default():

    bgc = MibigBGC(1, "BGC0000001", "BGC0000001", "Polyketide")
    assert bgc.id == 1
    assert bgc.strain == "BGC0000001"
    assert bgc.name == "BGC0000001"
    assert bgc.product_prediction == "Polyketide"
