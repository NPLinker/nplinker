from __future__ import annotations
import json
import pytest
from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.genomics import add_bgc_to_gcf
from nplinker.genomics import add_strain_to_bgc
from nplinker.genomics import generate_mappings_genome_id_bgc_id
from nplinker.genomics import get_mibig_from_gcf
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
def bgcs() -> list[BGC]:
    """Return a list of BGC objects."""
    return [BGC("BGC_01", "NPR"), BGC("BGC_02", "Alkaloid"), BGC("BGC_03", "Polyketide")]


def test_add_strain_to_bgc(bgcs):
    """Test add_strain_to_bgc function."""
    strain1 = Strain("STRAIN_01")
    strain1.add_alias("BGC_01")
    strain2 = Strain("STRAIN_02")
    strain2.add_alias("BGC_02")
    strain2.add_alias("BGC_02_1")
    strain3 = Strain("STRAIN_03")
    strains = StrainCollection()
    strains.add(strain1)
    strains.add(strain2)
    strains.add(strain3)

    bgc_with_strain, bgc_without_strain = add_strain_to_bgc(strains, bgcs)

    assert len(bgc_with_strain) == 2
    assert len(bgc_without_strain) == 1
    assert bgc_with_strain == [bgcs[0], bgcs[1]]
    assert bgc_without_strain == [bgcs[2]]
    assert bgc_with_strain[0].strain == strain1
    assert bgc_with_strain[1].strain == strain2
    assert bgc_without_strain[0].strain is None


def test_add_strain_to_bgc_error(bgcs):
    """Test add_strain_to_bgc function error."""
    strain1 = Strain("STRAIN_01")
    strain1.add_alias("BGC_01")
    strain2 = Strain("STRAIN_02")
    strain2.add_alias("BGC_01")
    strains = StrainCollection()
    strains.add(strain1)
    strains.add(strain2)

    with pytest.raises(ValueError) as e:
        add_strain_to_bgc(strains, bgcs)

    assert "Multiple strain objects found for BGC id 'BGC_01'" in e.value.args[0]


def test_add_bgc_to_gcf(bgcs):
    """Test add_bgc_to_gcf function."""
    gcf1 = GCF("1")
    gcf1.bgc_ids = {"BGC_01", "BGC_02"}
    gcf2 = GCF("2")
    gcf2.bgc_ids = {"BGC_03", "BGC_missing_01"}
    gcf3 = GCF("3")
    gcf3.bgc_ids = {"BGC_missing_02", "BGC_missing_03"}
    gcfs = [gcf1, gcf2, gcf3]

    gcf_with_bgc, gcf_without_bgc, gcf_missing_bgc = add_bgc_to_gcf(bgcs, gcfs)

    assert len(gcf_with_bgc) == 2
    assert len(gcf_without_bgc) == 1
    assert len(gcf_missing_bgc) == 2
    assert gcf_with_bgc == [gcf1, gcf2]
    assert gcf_without_bgc == [gcf3]
    assert gcf_missing_bgc == {gcf2: {"BGC_missing_01"}, gcf3: {"BGC_missing_02", "BGC_missing_03"}}
    assert gcf_with_bgc[0].bgcs == {bgcs[0], bgcs[1]}
    assert gcf_with_bgc[1].bgcs == {bgcs[2]}
    assert gcf_without_bgc[0].bgcs == set()


def test_get_mibig_from_gcf():
    """Test get_mibig_from_gcf function."""
    bgc1 = BGC("BGC_01", "NPR")
    bgc1.strain = Strain("BGC_01")
    bgc2 = BGC("BGC_02", "Alkaloid")
    bgc2.strain = Strain("BGC_02")
    bgc3 = BGC("antismash_c", "Polyketide")
    bgc3.strain = Strain("strain_01")
    gcf1 = GCF("1")
    gcf1.add_bgc(bgc1)
    gcf2 = GCF("2")
    gcf2.add_bgc(bgc2)
    gcf2.add_bgc(bgc3)
    gcfs = [gcf1, gcf2]

    mibig_bgcs_in_use, mibig_strains_in_use = get_mibig_from_gcf(gcfs)

    assert len(mibig_bgcs_in_use) == 2
    assert len(mibig_strains_in_use) == 2
    assert bgc3 not in mibig_bgcs_in_use
    assert bgc3.strain not in mibig_strains_in_use
