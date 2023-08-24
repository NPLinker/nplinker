import pytest
from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain


@pytest.fixture()
def bgc_with_strain():
    bgc = BGC("S0001", "NPR")
    bgc.strain = Strain("strain001")
    yield bgc


@pytest.fixture()
def bgc_without_strain():
    bgc = BGC("S002", "NPR")
    yield bgc


def test_default():
    gcf = GCF("1")
    assert gcf.gcf_id == "1"
    assert gcf.bgcs == set()
    assert isinstance(gcf.strains, StrainCollection)
    assert len(gcf.strains) == 0
    assert gcf.bigscape_class is None


def test_add_bgc(bgc_with_strain):
    gcf = GCF("1")
    gcf.add_bgc(bgc_with_strain)
    assert gcf.bgcs == {bgc_with_strain}
    assert bgc_with_strain.parents == {gcf}


def test_detach_bgc(bgc_with_strain):
    gcf = GCF("1")
    gcf.add_bgc(bgc_with_strain)
    gcf.detach_bgc(bgc_with_strain)
    assert gcf.bgcs == set()
    assert len(gcf.strains) == 0
    assert bgc_with_strain.parents == set()


def test_add_bgc_wo_strain(bgc_without_strain, caplog):
    gcf = GCF("1")
    gcf.add_bgc(bgc_without_strain)
    assert gcf.gcf_id == "1"
    assert gcf.bgcs == {bgc_without_strain}
    assert len(gcf.strains) == 0
    assert "No strain specified for the BGC" in caplog.text


def test_has_strain(bgc_with_strain):
    gcf = GCF("1")
    gcf.add_bgc(bgc_with_strain)
    assert gcf.has_strain(Strain("strain001")) is True
    assert gcf.has_strain(Strain("strain002")) is False

def test_has_mibig_only():
    mibig_bgc = BGC("BGC0000001", "NPR")
    nonmibig_bgc = BGC("S0001", "NPR")
    gcf = GCF("1")
    gcf.add_bgc(mibig_bgc)
    assert gcf.has_mibig_only() is True
    gcf.detach_bgc(mibig_bgc)
    gcf.add_bgc(nonmibig_bgc)
    assert gcf.has_mibig_only() is False
    gcf.add_bgc(mibig_bgc)
    assert gcf.has_mibig_only() is False
