import pytest
from nplinker.genomics import GCF
from nplinker.genomics.abc import GCFLoaderBase
from nplinker.genomics.bigscape import BigscapeV2GCFLoader
from .. import DATA_DIR


class TestBigscapelGCFLoader:
    @pytest.fixture
    def loader(self):
        db_file = DATA_DIR / "bigscape" / "mix" / "data_sqlite.db"
        loader = BigscapeV2GCFLoader(db_file)
        yield loader

    def test_abc(self, loader):
        assert issubclass(BigscapeV2GCFLoader, GCFLoaderBase)
        assert isinstance(loader, GCFLoaderBase)

    def test_init(self, loader):
        assert loader.db_file == str(
            DATA_DIR / "bigscape" / "mix" / "data_sqlite.db"
        )

    @pytest.mark.parametrize(
        "keep_mibig_only, keep_singleton, expected",
        [(False, False, 1), (True, False, 2), (False, True, 2), (True, True, 4)],
    )
    def test_get_gcfs(self, loader, keep_mibig_only, keep_singleton, expected):
        gcfs = loader.get_gcfs(keep_mibig_only, keep_singleton)
        assert isinstance(gcfs, list)
        assert len(gcfs) == expected
        assert isinstance(gcfs[0], GCF)

    def test_parse_gcf(self, loader):
        gcf_list = BigscapeV2GCFLoader._parse_gcf(loader.db_file)
        assert isinstance(gcf_list, list)
        assert len(gcf_list) == 4
        for gcf in gcf_list:
            assert isinstance(gcf, GCF)
