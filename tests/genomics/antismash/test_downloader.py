import pytest
from nplinker.genomics.antismash import download_and_extract_antismash_metadata

class TestDownloadAndExtractAntismashData():

    def test_default(self, tmp_path):
        download_root = tmp_path / "download"
        extract_root = tmp_path / "extracted"
        download_root.mkdir()
        extract_root.mkdir()
        refseq_assembly_id = "GCF_004339725.1"
        download_and_extract_antismash_metadata(refseq_assembly_id, download_root, extract_root)
        archive = download_root / "GCF_004339725.1.zip"
        extracted_folder = extract_root / "antismash/GCF_004339725.1"
        assert archive.exists()
        assert archive.is_file()
        assert extracted_folder.exists()

    def test_error_same_path(self, tmp_path):
        refseq_assembly_id = "GCF_004339725.1"
        with pytest.raises(ValueError) as e:
            download_and_extract_antismash_metadata(refseq_assembly_id, tmp_path, tmp_path)
        assert e.value.args[
            0] == "Identical path of download directory and extract directory"

    def test_error_nonempty_path(self, tmp_path):
        refseq_assembly_id = "GCF_004339725.1"
        nonempty_path = tmp_path / f"extracted/antismash/{refseq_assembly_id}/subdir"
        nonempty_path.mkdir(parents=True)
        with pytest.raises(ValueError) as e:
            download_and_extract_antismash_metadata(
                refseq_assembly_id,
                tmp_path,
                tmp_path / "extracted")
        assert "Nonempty directory" in e.value.args[0]
