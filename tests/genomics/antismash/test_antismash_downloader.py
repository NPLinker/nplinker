import pytest
import urllib
from nplinker.genomics.antismash import download_and_extract_antismash_data
from nplinker.utils import list_files, extract_archive


class TestDownloadAndExtractAntismashData():

    antismash_id = "GCF_004339725.1"

    def test_default(self, tmp_path):
        download_root = tmp_path / "download"
        download_root.mkdir()
        extract_root = tmp_path / "extracted"
        extract_root.mkdir()
        original_extract_root = tmp_path / "original"
        original_extract_root.mkdir()
        download_and_extract_antismash_data(self.antismash_id, download_root,
                                            extract_root)
        archive = download_root / "GCF_004339725.1.zip"
        extracted_folder = extract_root / "antismash" / "GCF_004339725.1"
        extracted_files = list_files(extracted_folder, keep_parent=False)
        # extract zip folder without removing any files
        extract_archive(archive, original_extract_root)
        original_expected_files = list_files(original_extract_root,
                                             suffix=(".json", ".gbk"),
                                             keep_parent=False)
        assert archive.exists()
        assert archive.is_file()
        assert extracted_folder.exists()
        assert extracted_files.sort() == original_expected_files.sort()

    def test_error_same_path(self, tmp_path):
        with pytest.raises(ValueError) as e:
            download_and_extract_antismash_data(self.antismash_id, tmp_path,
                                                tmp_path)
        assert e.value.args[
            0] == "Identical path of download directory and extract directory"

    def test_error_nonempty_path(self, tmp_path):
        nonempty_path = tmp_path / "extracted" / "antismash" / f"{self.antismash_id}" / "subdir"
        nonempty_path.mkdir(parents=True)
        with pytest.raises(ValueError) as e:
            download_and_extract_antismash_data(self.antismash_id, tmp_path,
                                                tmp_path / "extracted")
        assert "Nonempty directory" in e.value.args[0]

    def test_broken_url(self, tmp_path):
        broken_antismash_id = 'broken_id'
        download_root = tmp_path / "download"
        download_root.mkdir()
        extract_root = tmp_path / "extracted"
        extract_root.mkdir()
        with pytest.raises(urllib.error.HTTPError):
            download_and_extract_antismash_data(broken_antismash_id, download_root,
                                                extract_root)
            archive = download_root / broken_antismash_id / ".zip"
            extracted_folder = extract_root / "antismash" / broken_antismash_id
            assert not archive.exists()
            assert not archive.is_file()
            assert not extracted_folder.exists()
