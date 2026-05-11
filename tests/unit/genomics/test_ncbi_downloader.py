from unittest.mock import patch
import httpx
import pytest
from nplinker.genomics.antismash.ncbi_downloader import download_and_extract_ncbi_genome


@pytest.fixture
def download_root(tmp_path):
    return tmp_path / "download"


@pytest.fixture
def extract_root(tmp_path):
    return tmp_path / "extracted"


def test_download_and_extract_ncbi_genome_success(download_root, extract_root):
    assembly_accession = "GCF_000514775.1"

    genome_path = download_and_extract_ncbi_genome(assembly_accession, download_root, extract_root)

    assert genome_path == extract_root / "ncbi_genomes" / f"{assembly_accession}.gbff"
    assert not (extract_root / "ncbi_genomes" / "md5sum.txt").exists()
    assert not (extract_root / "ncbi_genomes" / "README.md").exists()
    assert not (extract_root / "ncbi_genomes" / "ncbi_dataset").exists()


def test_download_and_extract_ncbi_genome_max_retries(download_root, extract_root):
    assembly_accession = "GCF_000514775.1"

    with patch(
        "nplinker.genomics.antismash.ncbi_downloader.download_url",
        side_effect=httpx.ReadTimeout("Download failed"),
    ):
        with pytest.raises(httpx.ReadTimeout, match="Maximum download retries"):
            download_and_extract_ncbi_genome(
                assembly_accession, download_root, extract_root, max_attempts=1
            )


def test_download_and_extract_ncbi_genome_invalid_accession(download_root, extract_root):
    assembly_accession = "invalid_ref_seq_id"

    with pytest.raises(ValueError, match="Not a valid genome assembly accession"):
        download_and_extract_ncbi_genome(assembly_accession, download_root, extract_root)
