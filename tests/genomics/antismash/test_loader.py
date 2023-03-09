import pytest
from nplinker.genomics import BGC
from nplinker.genomics import BGCLoaderBase
from nplinker.genomics.antismash import AntismashBGCLoader
from nplinker.genomics.antismash import parse_bgc_genbank
from ... import DATA_DIR


class TestAntismashBGCLoader:

    @pytest.fixture
    def loader(self):
        data_dir = str(DATA_DIR / "antismash")
        loader = AntismashBGCLoader(data_dir)
        yield loader

    def test_abc(self, loader):
        assert issubclass(AntismashBGCLoader, BGCLoaderBase)
        assert isinstance(loader, BGCLoaderBase)

    def test_init(self, loader):
        assert loader.data_dir == str(DATA_DIR / "antismash")

    def test_get_files(self, loader):
        bgc_files = loader.get_files()
        assert isinstance(bgc_files, dict)
        assert len(bgc_files) == 44
        assert "NZ_AZWB01000005.region001" in bgc_files
        assert "NZ_AZWS01000001.region001" in bgc_files
        assert "GCF_000514855.1" not in bgc_files
        assert "GCF_000514515.1" not in bgc_files

    def test_parse_data_dir(self):
        data_dir = DATA_DIR / "antismash"
        bgc_files = AntismashBGCLoader._parse_data_dir(str(data_dir))
        assert isinstance(bgc_files, dict)
        assert len(bgc_files) == 44
        assert "NZ_AZWB01000005.region001" in bgc_files
        assert "NZ_AZWS01000001.region001" in bgc_files
        assert "GCF_000514855.1" not in bgc_files
        assert "GCF_000514515.1" not in bgc_files
        assert bgc_files["NZ_AZWB01000005.region001"] == str(
            data_dir / "GCF_000514515.1" / "NZ_AZWB01000005.region001.gbk")

    def test_get_bgcs(self, loader):
        bgcs = loader.get_bgcs()
        assert isinstance(bgcs, dict)
        assert len(bgcs) == 44
        assert isinstance(bgcs["NZ_AZWB01000005.region001"], BGC)
        assert isinstance(bgcs["NZ_AZWS01000001.region001"], BGC)
        assert bgcs.get("GCF_000514855.1", "NotExist") == "NotExist"
        assert bgcs.get("GCF_000514515.1", "NotExist") == "NotExist"


def test_parse_bgc_genbank():
    gbk_file = str(DATA_DIR / "antismash" / "GCF_000514515.1" /
                   "NZ_AZWB01000005.region001.gbk")
    bgc = parse_bgc_genbank(gbk_file)
    assert isinstance(bgc, BGC)
    assert bgc.bgc_id == "NZ_AZWB01000005.region001"
    assert bgc.product_prediction == ["NRPS", "lanthipeptide"]
    assert "Salinispora pacifica CNT029 B170DRAFT_scaffold" in bgc.description
    assert bgc.antismash_id == "NZ_AZWB01000005"
    assert bgc.antismash_file == gbk_file
    assert bgc.antismash_region == "1"
    assert bgc.smiles == [
        "NC([*])C(=O)NC([*])C(=O)NC(CO)C(=O)NC(Cc1ccccc1)C(=O)NCC(=O)O"
    ]


def test_parse_bgc_genbank_error():
    gbk_file = str(DATA_DIR / "fake_antismash.region001.gbk")
    with pytest.raises(ValueError) as e:
        parse_bgc_genbank(gbk_file)
    assert "Not found product prediction in antiSMASH Genbank file" in e.value.args[
        0]
