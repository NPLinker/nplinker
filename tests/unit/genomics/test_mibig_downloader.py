import pytest
from nplinker.genomics import mibig


class TestDownloadAndExtractMibigMetadata:
    @pytest.mark.parametrize(
        "version, expected",
        [
            ["1.4", "mibig_json_1.4.tar.gz"],
            ["3.1", "mibig_json_3.1.tar.gz"],
            ["4.0", "mibig_json_4.0.tar.gz"],
        ],
    )
    def test_version(self, tmp_path, version, expected):
        download_path = tmp_path / "download"
        extract_path = tmp_path / "metadata"
        download_path.mkdir()
        extract_path.mkdir()
        mibig.download_and_extract_mibig_metadata(download_path, extract_path, version=version)
        archive = download_path / expected
        metadata = extract_path / "BGC0000002.json"
        assert archive.exists()
        assert archive.is_file()
        assert metadata.exists()
        assert metadata.is_file()

    def test_error_same_path(self, tmp_path):
        with pytest.raises(
            ValueError, match="Identical path of download directory and extract directory"
        ):
            mibig.download_and_extract_mibig_metadata(tmp_path, tmp_path)

    def test_error_nonempty_path(self, tmp_path):
        nonempty_path = tmp_path / "metadata" / "subdir"
        nonempty_path.mkdir(parents=True)

        with pytest.raises(ValueError, match="Nonempty directory"):
            mibig.download_and_extract_mibig_metadata(tmp_path, nonempty_path.parent)
