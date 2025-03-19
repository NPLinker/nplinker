import logging
import os
import shutil
import time
from os import PathLike
from pathlib import Path
from typing import Optional
from nplinker.utils import check_md5
from nplinker.utils import download_url
from nplinker.utils import extract_archive


logger = logging.getLogger(__name__)


def download_and_extract_ncbi_genome(
    genome_assembly_acc: str,
    download_root: str | PathLike,
    extract_root: str | PathLike,
    max_attempts: int = 10,
) -> Optional[Path]:
    """Downloads and extracts an NCBI dataset for a given genome RefSeq ID.

    This function retrieves a dataset from the NCBI database using the provided
    RefSeq ID. It retries the download process up to a specified maximum number
    of attempts in case of errors. The function verifies the integrity of the
    downloaded files using MD5 checksums, extracts the dataset, and renames the
    GenBank file for easier access. Unnecessary files are removed after successful
    processing.

    Args:
        genome_assembly_acc (str): The NCBI accession of the genome assembly to be downloaded.
        download_root (str | PathLike): The directory where the dataset will be downloaded.
        extract_root (str | PathLike): The directory where the dataset will be extracted.
        max_attempts (int): The maximum number of download attempts. Defaults to 10.

    Returns:
        Optional[Path]: The path to the extracted GenBank file if successful, otherwise None.

    Raises:
        RuntimeError: If the maximum number of retries is reached and the dataset
            could not be successfully downloaded and extracted.
    """
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{genome_assembly_acc}/download?include_annotation_type=GENOME_GB"
    )

    download_root = Path(download_root)
    extract_path = Path(extract_root) / "ncbi_genomes"
    filename = f"ncbi_{genome_assembly_acc}.zip"

    extract_path.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, max_attempts + 1):
        try:
            download_url(url, download_root, filename)
            archive = download_root / filename
            break
        except Exception as e:
            logger.warning(f"Attempt {attempt}/{max_attempts} failed to download {url}. Error: {e}")
            if attempt < max_attempts:
                time.sleep(2)
    else:
        raise RuntimeError(
            f"Maximum download retries ({max_attempts}) reached for {url}. Download failed."
        )

    extract_archive(archive, extract_path)
    verify_ncbi_dataset_md5_sums(extract_path)

    # Move and rename GenBank file
    genbank_path = extract_path / "ncbi_dataset" / "data" / genome_assembly_acc / "genomic.gbff"
    new_genbank_path = extract_path / f"{genome_assembly_acc}.gbff"
    genbank_path.rename(new_genbank_path)

    # Delete unnecessary files
    shutil.rmtree(extract_path / "ncbi_dataset")
    os.remove(extract_path / "md5sum.txt")
    os.remove(extract_path / "README.md")

    return new_genbank_path


def verify_ncbi_dataset_md5_sums(extract_path: PathLike) -> bool:
    """Verify the integrity of files in a specified directory using MD5 checksums.

    This function reads an "md5sum.txt" file located in the given extraction path,
    which contains MD5 checksums and corresponding file names. It then computes
    the MD5 checksum for each file and compares it with the expected value. If any
    file's checksum does not match, a `ValueError` is raised.

    Args:
        extract_path (PathLike): Path to the directory containing the files and
            the "md5sum.txt" file.

    Returns:
        bool: True if all files pass the MD5 checksum verification.

    Raises:
        ValueError: If the MD5 checksum of any file does not match the expected value.
    """
    with open(extract_path / "md5sum.txt", "r") as f:
        for line in f:
            md5sum, file_name = line.strip().split()
            file_path = extract_path / file_name
            if not check_md5(file_path, md5sum):
                raise ValueError(f"MD5 checksum mismatch for {file_path}")
