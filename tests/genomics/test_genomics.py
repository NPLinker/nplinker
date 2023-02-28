from __future__ import annotations
import pytest
from nplinker.genomics import BGC
from nplinker.genomics import filter_mibig_only_gcf
from nplinker.genomics import GCF
from nplinker.genomics import get_bgcs_from_gcfs
from nplinker.genomics import get_strains_from_bgcs
from nplinker.genomics import map_bgc_to_gcf
from nplinker.genomics import map_strain_to_bgc
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain


@pytest.fixture
def strain_collection() -> StrainCollection:
    sc = StrainCollection()
    sc.add(Strain("SAMPLE0001"))
    sc.add(Strain("BGC0000001"))

    strain = Strain("EXAMPLE002")
    strain.add_alias("SAMPLE0002")
    sc.add(strain)

    strain = Strain("EXAMPLE003")
    strain.add_alias("BGC0000002")
    sc.add(strain)
    return sc


@pytest.fixture
def bgc_list() -> list[BGC]:
    return [
        BGC("SAMPLE0001", ["NPR"]),
        BGC("SAMPLE0002", ["Alkaloid"]),
        BGC("BGC0000001", ["Polyketide"]),
        BGC("BGC0000002", ["Terpene"])
    ]


@pytest.fixture
def bgc_list_error() -> list[BGC]:
    return [
        BGC("SAMPLE0003", ["NPR"]),
        BGC("BGC0000003", ["Polyketide"]),
    ]


@pytest.fixture
def gcf_list() -> list[GCF]:
    gcf1 = GCF("1")
    gcf1.bgc_ids |= set(("SAMPLE0001", "SAMPLE0002"))
    gcf2 = GCF("2")
    gcf2.bgc_ids |= set(("BGC0000001", "BGC0000002"))
    return [gcf1, gcf2]


@pytest.fixture
def gcf_list_error() -> list[GCF]:
    gcf1 = GCF("1")
    gcf1.bgc_ids |= set(("SAMPLE0001", "SAMPLE0003"))
    return [gcf1]


def test_map_strain_to_bgc(strain_collection, bgc_list):
    for bgc in bgc_list:
        assert bgc.strain is None
    map_strain_to_bgc(strain_collection, bgc_list)
    for bgc in bgc_list:
        assert bgc.strain is not None
    assert bgc_list[0].strain.id == "SAMPLE0001"
    assert bgc_list[1].strain.id == "EXAMPLE002"
    assert bgc_list[2].strain.id == "BGC0000001"
    assert bgc_list[3].strain.id == "EXAMPLE003"


def test_map_strain_to_bgc_error(strain_collection, bgc_list_error):
    for bgc in bgc_list_error:
        assert bgc.strain is None
    with pytest.raises(KeyError) as e:
        map_strain_to_bgc(strain_collection, bgc_list_error)
    assert "Strain id SAMPLE0003 from BGC object SAMPLE0003 not found" in e.value.args[
        0]


def test_map_bgc_to_gcf(bgc_list, gcf_list):
    assert gcf_list[0].bgc_ids == set(("SAMPLE0001", "SAMPLE0002"))
    assert gcf_list[1].bgc_ids == set(("BGC0000001", "BGC0000002"))
    assert len(gcf_list[0].bgcs) == 0
    assert len(gcf_list[1].bgcs) == 0
    map_bgc_to_gcf(bgc_list, gcf_list)
    assert gcf_list[0].bgc_ids == set(("SAMPLE0001", "SAMPLE0002"))
    assert gcf_list[1].bgc_ids == set(("BGC0000001", "BGC0000002"))
    assert len(gcf_list[0].bgcs) == 2
    assert len(gcf_list[1].bgcs) == 2
    assert gcf_list[0].bgcs == set(bgc_list[:2])
    assert gcf_list[1].bgcs == set(bgc_list[2:])


def test_map_bgc_to_gcf_error(bgc_list, gcf_list_error):
    assert gcf_list_error[0].bgc_ids == set(("SAMPLE0001", "SAMPLE0003"))
    assert len(gcf_list_error[0].bgcs) == 0
    with pytest.raises(KeyError) as e:
        map_bgc_to_gcf(bgc_list, gcf_list_error)
    assert "BGC id SAMPLE0003 from GCF object 1 not found" in e.value.args[0]


def test_filter_mibig_only_gcf(bgc_list, gcf_list):
    map_bgc_to_gcf(bgc_list, gcf_list)
    gcfs = filter_mibig_only_gcf(gcf_list)
    assert len(gcfs) == 1
    assert gcfs[0].gcf_id == "1"


def test_get_bgcs_from_gcfs(bgc_list, gcf_list):
    map_bgc_to_gcf(bgc_list, gcf_list)
    bgcs = get_bgcs_from_gcfs(gcf_list)
    assert isinstance(bgcs, list)
    assert len(bgcs) == 4
    for i in bgcs:
        assert isinstance(i, BGC)


def test_get_strains_from_bgcs(strain_collection, bgc_list):
    map_strain_to_bgc(strain_collection, bgc_list)
    strains = get_strains_from_bgcs(bgc_list)
    assert isinstance(strains, StrainCollection)
    for strain in strain_collection:
        assert strain in strains
    # assert strains == strain_collection # use it when issue #113 is solved
