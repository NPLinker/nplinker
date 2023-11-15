import pytest
from nplinker.genomics import GCF
from nplinker.genomics.abc import GCFLoaderBase
from nplinker.genomics.bigscape import BigscapeGCFLoader
from .. import DATA_DIR


class TestBigscapelGCFLoader:
    @pytest.fixture
    def loader(self):
        cluster_file = DATA_DIR / "bigscape" / "mix" / "mix_clustering_c0.30.tsv"
        loader = BigscapeGCFLoader(cluster_file)
        yield loader

    def test_abc(self, loader):
        assert issubclass(BigscapeGCFLoader, GCFLoaderBase)
        assert isinstance(loader, GCFLoaderBase)

    def test_init(self, loader):
        assert loader.cluster_file == str(
            DATA_DIR / "bigscape" / "mix" / "mix_clustering_c0.30.tsv"
        )

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
        gcf_list = BigscapeGCFLoader._parse_gcf(loader.cluster_file)  # noqa
        assert isinstance(gcf_list, list)
        assert len(gcf_list) == 5
        for gcf in gcf_list:
            assert isinstance(gcf, GCF)
