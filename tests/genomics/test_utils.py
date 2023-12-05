from __future__ import annotations
import json
import pytest
from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.genomics import add_bgc_to_gcf
from nplinker.genomics import add_strain_to_bgc
from nplinker.genomics import generate_mappings_genome_id_bgc_id
from nplinker.genomics import get_bgcs_from_gcfs
from nplinker.genomics import get_strains_from_bgcs
from nplinker.globals import GENOME_BGC_MAPPINGS_FILENAME
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection
from .. import DATA_DIR


def test_generate_mappings_genome_id_bgc_id(tmp_path):
    """Test generate_mappings_genome_id_bgc_id function."""
    bgc_dir = DATA_DIR / "antismash"

    # using default output file path
    generate_mappings_genome_id_bgc_id(bgc_dir)
    # using custom output file path
    generate_mappings_genome_id_bgc_id(bgc_dir, tmp_path / GENOME_BGC_MAPPINGS_FILENAME)

    # read both files
    with open(bgc_dir / GENOME_BGC_MAPPINGS_FILENAME) as f:
        mappings = json.load(f)
    with open(tmp_path / GENOME_BGC_MAPPINGS_FILENAME) as f:
        mappings_with_outfile = json.load(f)

    # check if both files are the same
    assert mappings == mappings_with_outfile

    # then check the content
    assert len(mappings["mappings"]) == 2

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


def test_generate_mappings_genome_id_bgc_id_empty_dir(tmp_path, caplog):
    """Test generate_mappings_genome_id_bgc_id function with empty dir."""
    # prepare dir and file
    bgc_dir = tmp_path / "GCF_1"
    bgc_file = bgc_dir / "BGC_1.gbk"
    bgc_dir.mkdir()
    bgc_file.touch()
    empty_dir = tmp_path / "empty_dir"
    empty_dir.mkdir()

    generate_mappings_genome_id_bgc_id(tmp_path)
    assert "No BGC files found" in caplog.text

    with open(tmp_path / GENOME_BGC_MAPPINGS_FILENAME) as f:
        mappings = json.load(f)
    assert len(mappings["mappings"]) == 1
    assert mappings["mappings"][0]["genome_ID"] == "GCF_1"
    assert mappings["mappings"][0]["BGC_ID"] == ["BGC_1"]


@pytest.fixture
def strain_collection() -> StrainCollection:
    """Return a StrainCollection object."""
    sc = StrainCollection()

    strain = Strain("STRAIN_01")
    strain.add_alias("BGC_01")
    sc.add(strain)

    strain = Strain("STRAIN_02")
    strain.add_alias("BGC_02")
    strain.add_alias("BGC_02_1")
    sc.add(strain)

    strain = Strain("SAMPLE_BGC_03")
    sc.add(strain)

    return sc


@pytest.fixture
def bgc_list() -> list[BGC]:
    """Return a list of BGC objects."""
    return [BGC("BGC_01", "NPR"), BGC("BGC_02", "Alkaloid"), BGC("SAMPLE_BGC_03", "Polyketide")]


@pytest.fixture
def gcf_list() -> list[GCF]:
    """Return a list of GCF objects."""
    gcf1 = GCF("1")
    gcf1.bgc_ids |= {"BGC_01"}
    gcf2 = GCF("2")
    gcf2.bgc_ids |= {"BGC_02", "SAMPLE_BGC_03"}
    return [gcf1, gcf2]


@pytest.fixture
def gcf_list_error() -> list[GCF]:
    """Return a list of GCF objects for testing errors."""
    gcf1 = GCF("1")
    gcf1.bgc_ids |= {"SAMPLE_BGC_03", "BGC_04"}
    return [gcf1]


def test_add_strain_to_bgc(strain_collection, bgc_list):
    """Test add_strain_to_bgc function."""
    for bgc in bgc_list:
        assert bgc.strain is None
    add_strain_to_bgc(strain_collection, bgc_list)
    for bgc in bgc_list:
        assert bgc.strain is not None
    assert bgc_list[0].strain.id == "STRAIN_01"
    assert bgc_list[1].strain.id == "STRAIN_02"
    assert bgc_list[2].strain.id == "SAMPLE_BGC_03"


def test_add_strain_to_bgc_error(strain_collection):
    """Test add_strain_to_bgc function error."""
    bgcs = [BGC("BGC_04", "NPR")]
    with pytest.raises(ValueError) as e:
        add_strain_to_bgc(strain_collection, bgcs)
    assert "Strain id 'BGC_04' from BGC object 'BGC_04' not found" in e.value.args[0]


def test_add_bgc_to_gcf(bgc_list, gcf_list):
    """Test add_bgc_to_gcf function."""
    assert gcf_list[0].bgc_ids == {"BGC_01"}
    assert gcf_list[1].bgc_ids == {"BGC_02", "SAMPLE_BGC_03"}
    assert len(gcf_list[0].bgcs) == 0
    assert len(gcf_list[1].bgcs) == 0
    add_bgc_to_gcf(bgc_list, gcf_list)
    assert gcf_list[0].bgc_ids == {"BGC_01"}
    assert gcf_list[1].bgc_ids == {"BGC_02", "SAMPLE_BGC_03"}
    assert len(gcf_list[0].bgcs) == 1
    assert len(gcf_list[1].bgcs) == 2
    assert gcf_list[0].bgcs == set(bgc_list[:1])
    assert gcf_list[1].bgcs == set(bgc_list[1:])


def test_add_bgc_to_gcf_error(bgc_list, gcf_list_error):
    """Test add_bgc_to_gcf function error."""
    assert gcf_list_error[0].bgc_ids == {"SAMPLE_BGC_03", "BGC_04"}
    assert len(gcf_list_error[0].bgcs) == 0
    with pytest.raises(KeyError) as e:
        add_bgc_to_gcf(bgc_list, gcf_list_error)
    assert "BGC id 'BGC_04' from GCF object '1' not found" in e.value.args[0]


def test_get_bgcs_from_gcfs(bgc_list, gcf_list):
    """Test get_bgcs_from_gcfs function."""
    add_bgc_to_gcf(bgc_list, gcf_list)
    bgcs = get_bgcs_from_gcfs(gcf_list)
    assert isinstance(bgcs, list)
    assert len(bgcs) == 3
    for i in bgcs:
        assert isinstance(i, BGC)


def test_get_strains_from_bgcs(strain_collection, bgc_list):
    """Test get_strains_from_bgcs function."""
    add_strain_to_bgc(strain_collection, bgc_list)
    strains = get_strains_from_bgcs(bgc_list)
    assert isinstance(strains, StrainCollection)
    assert strains == strain_collection
