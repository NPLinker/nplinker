import pytest
from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.strain import Strain
from nplinker.strain import StrainCollection


@pytest.fixture()
def bgc_with_strain():
    """Return a BGC object that has a Strain."""
    bgc = BGC("S0001", "NPR")
    bgc.strain = Strain("strain001")
    yield bgc


@pytest.fixture()
def bgc_without_strain():
    """Return a BGC object that does not have Strain."""
    bgc = BGC("S002", "NPR")
    yield bgc


def test_init():
    """Test the initialization of GCF."""
    gcf = GCF("1")
    assert gcf.id == "1"
    assert gcf.bgcs == set()
    assert isinstance(gcf.strains, StrainCollection)
    assert len(gcf.strains) == 0
    assert gcf.bigscape_class is None


def test_str_repr_():
    """Test __str__ and __repr__ method."""
    gcf = GCF("1")
    assert str(gcf) == "GCF(id=1, #BGC_objects=0, #bgc_ids=0,#strains=0)."
    assert repr(gcf) == "GCF(id=1, #BGC_objects=0, #bgc_ids=0,#strains=0)."


def test_eq():
    """Test __eq__ method."""
    gcf1 = GCF("1")
    gcf2 = GCF("1")
    assert gcf1 == gcf2
    gcf3 = GCF("2")
    assert gcf1 != gcf3


def test_hash():
    """Test __hash__ method."""
    gcf1 = GCF("1")
    gcf2 = GCF("1")
    assert hash(gcf1) == hash(gcf2)
    gcf3 = GCF("2")
    assert hash(gcf1) != hash(gcf3)


def test_add_bgc(bgc_with_strain):
    """Test add_bgc method."""
    gcf = GCF("1")
    gcf.add_bgc(bgc_with_strain)
    assert gcf.bgcs == {bgc_with_strain}
    assert bgc_with_strain.parents == {gcf}


def test_detach_bgc(bgc_with_strain):
    """Test detach_bgc method."""
    gcf = GCF("1")
    gcf.add_bgc(bgc_with_strain)
    gcf.detach_bgc(bgc_with_strain)
    assert gcf.bgcs == set()
    assert len(gcf.strains) == 0
    assert bgc_with_strain.parents == set()


def test_add_bgc_wo_strain(bgc_without_strain, caplog):
    """Test add_bgc method with a BGC that does have strain."""
    gcf = GCF("1")
    gcf.add_bgc(bgc_without_strain)
    assert gcf.id == "1"
    assert gcf.bgcs == {bgc_without_strain}
    assert len(gcf.strains) == 0
    assert "No strain specified for the BGC" in caplog.text


def test_has_strain(bgc_with_strain):
    """Test has_strain method."""
    gcf = GCF("1")
    gcf.add_bgc(bgc_with_strain)
    assert gcf.has_strain(Strain("strain001")) is True
    assert gcf.has_strain(Strain("strain002")) is False


def test_has_mibig_only():
    """Test has_mibig_only method."""
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


def test_is_singleton():
    """Test is_singleton method."""
    bgc1 = BGC("BGC0000001", "NPR")
    bgc2 = BGC("BGC0000002", "NPR")
    gcf = GCF("1")
    gcf.add_bgc(bgc1)
    assert gcf.is_singleton() is True
    gcf.add_bgc(bgc2)
    assert gcf.is_singleton() is False
    gcf.detach_bgc(bgc1)
    assert gcf.is_singleton() is True
