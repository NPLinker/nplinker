from __future__ import annotations
import json
import pytest
from nplinker.genomics import BGC
from nplinker.genomics import filter_mibig_only_gcf
from nplinker.genomics import GCF
from nplinker.genomics import generate_genome_bgc_mappings_file
from nplinker.genomics import GENOME_BGC_MAPPINGS_FILENAME
from nplinker.genomics import get_bgcs_from_gcfs
from nplinker.genomics import get_strains_from_bgcs
from nplinker.genomics import map_bgc_to_gcf
from nplinker.genomics import map_strain_to_bgc
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain
from .. import DATA_DIR


def test_generate_genome_bgc_mappings_file():
    bgc_dir = DATA_DIR / "antismash"

    generate_genome_bgc_mappings_file(bgc_dir)

    with open(bgc_dir / GENOME_BGC_MAPPINGS_FILENAME) as f:
        mappings = json.load(f)

    assert mappings["count"] == 2

    assert mappings["mappings"][0]["genome_ID"] == "GCF_000514515.1"
    assert len(mappings["mappings"][0]["BGC_ID"]) == 20
    for bgc_id in ["NZ_AZWB01000005.region001", "NZ_KI911412.1.region001"]:
        assert bgc_id in mappings["mappings"][0]["BGC_ID"]
    assert "GCF_000514515.1" not in mappings["mappings"][0]["BGC_ID"]

    assert mappings["mappings"][1]["genome_ID"] == "GCF_000514855.1"
    assert len(mappings["mappings"][1]["BGC_ID"]) == 24
    for bgc_id in ["NZ_AZWS01000001.region001", "NZ_KI911483.1.region001"]:
        assert bgc_id in mappings["mappings"][1]["BGC_ID"]
    assert "GCF_000514855.1" not in mappings["mappings"][1]["BGC_ID"]

    (bgc_dir / GENOME_BGC_MAPPINGS_FILENAME).unlink()


@pytest.fixture
def strain_collection() -> StrainCollection:
    sc = StrainCollection()
    sc.add(Strain("G-SAMPLE0001"))

    strain = Strain("S-EXAMPLE002")
    strain.add_alias("G-SAMPLE0002")
    sc.add(strain)

    sc.add(Strain("BGC0000001"))
    sc.add(Strain("BGC0000002"))
    return sc


@pytest.fixture
def bgc_genome_mapping() -> dict[str, str]:
    return {
        "SAMPLE0001": "G-SAMPLE0001",
        "SAMPLE0002": "G-SAMPLE0002",
        "BGC0000001": "BGC0000001",
        "BGC0000002": "BGC0000002"
    }


@pytest.fixture
def bgc_list() -> list[BGC]:
    return [
        BGC("SAMPLE0001", "NPR"),
        BGC("SAMPLE0002", "Alkaloid"),
        BGC("BGC0000001", "Polyketide"),
        BGC("BGC0000002", "Terpene")
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


def test_map_strain_to_bgc(strain_collection, bgc_list, bgc_genome_mapping):
    for bgc in bgc_list:
        assert bgc.strain is None
    map_strain_to_bgc(strain_collection, bgc_list, bgc_genome_mapping)
    for bgc in bgc_list:
        assert bgc.strain is not None
    assert bgc_list[0].strain.id == "G-SAMPLE0001"
    assert bgc_list[1].strain.id == "S-EXAMPLE002"
    assert bgc_list[2].strain.id == "BGC0000001"
    assert bgc_list[3].strain.id == "BGC0000002"


def test_map_strain_to_bgc_error(strain_collection):
    bgc_genome_mapping = {"SAMPLE0003": "G-SAMPLE0003"}
    bgcs = [BGC("BGC0000003", ["Polyketide"])]
    with pytest.raises(KeyError) as e:
        map_strain_to_bgc(strain_collection, bgcs, bgc_genome_mapping)
    assert "Not found BGC id BGC0000003 in BGC-genome mappings." in e.value.args[
        0]

    bgcs = [BGC("SAMPLE0003", ["NPR"])]
    with pytest.raises(KeyError) as e:
        map_strain_to_bgc(strain_collection, bgcs, bgc_genome_mapping)
    assert "Strain id G-SAMPLE0003 from BGC object SAMPLE0003 not found" in e.value.args[
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


def test_get_strains_from_bgcs(strain_collection, bgc_list,
                               bgc_genome_mapping):
    map_strain_to_bgc(strain_collection, bgc_list, bgc_genome_mapping)
    strains = get_strains_from_bgcs(bgc_list)
    assert isinstance(strains, StrainCollection)
    for strain in strain_collection:
        assert strain in strains
    # assert strains == strain_collection # use it when issue #113 is solved
