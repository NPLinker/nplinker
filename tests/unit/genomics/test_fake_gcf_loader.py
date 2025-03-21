import pytest
from nplinker.genomics import GCF
from nplinker.genomics.abc import GCFLoaderBase
from nplinker.genomics.fake_gcf_loader import FakeGCFLoader
from .. import DATA_DIR


class TestFakelGCFLoader:
    @pytest.fixture
    def loader(self):
        cluster_file = DATA_DIR / "fake_gcf.csv"
        loader = FakeGCFLoader(cluster_file)
        yield loader

    def test_abc(self, loader):
        assert issubclass(FakeGCFLoader, GCFLoaderBase)
        assert isinstance(loader, GCFLoaderBase)

    def test_init(self, loader):
        assert loader.cluster_file == str(DATA_DIR / "fake_gcf.csv")

    @pytest.mark.parametrize(
        "keep_mibig_only, keep_singleton, expected",
        [(False, False, 2), (True, False, 3), (False, True, 3), (True, True, 5)],
    )
    def test_get_gcfs(self, loader, keep_mibig_only, keep_singleton, expected):
        gcfs = loader.get_gcfs(keep_mibig_only, keep_singleton)
        assert isinstance(gcfs, list)
        assert len(gcfs) == expected
        assert isinstance(gcfs[0], GCF)

    def test_parse_gcf(self, loader):
        # It's not really needed to test private method
        gcf_list = FakeGCFLoader._parse_gcf(loader.cluster_file)  # noqa
        assert isinstance(gcf_list, list)
        assert len(gcf_list) == 5
        for gcf in gcf_list:
            assert isinstance(gcf, GCF)
