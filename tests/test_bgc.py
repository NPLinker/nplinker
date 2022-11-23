from nplinker.genomics import BGC


def test_default():

    bgc = BGC(1, "BGC0000001", ["Polyketide"])
    assert bgc.id == 1
    assert bgc.name == "BGC0000001"
    assert bgc.product_prediction == ["Polyketide"]
    assert bgc.is_mibig() is True
