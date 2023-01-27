from nplinker.genomics import BGC


def test_default():

    bgc = BGC("BGC0000001", ["Polyketide"])
    assert bgc.bgc_id == "BGC0000001"
    assert bgc.product_prediction == ["Polyketide"]
    assert bgc.is_mibig() is True
