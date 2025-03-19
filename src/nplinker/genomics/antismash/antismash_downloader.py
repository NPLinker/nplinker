from __future__ import annotations
import logging
import os
import shutil
from os import PathLike
from pathlib import Path
import requests
from nplinker.utils import download_and_extract_archive
from nplinker.utils import list_dirs
from nplinker.utils import list_files


logger = logging.getLogger(__name__)

# urls to be given to download antismash data
ANTISMASH_DB_DOWNLOAD_URL = "https://antismash-db.secondarymetabolites.org/output/{}/{}"
# The antiSMASH DBV2 is for the availability of the old version, better to keep it.
ANTISMASH_DBV2_DOWNLOAD_URL = "https://antismash-dbv2.secondarymetabolites.org/output/{}/{}"
# antismash api to download results from submitted jobs
ANTISMASH_API_DOWNLOAD_URL = "https://antismash.secondarymetabolites.org/upload/{}/{}"


def download_and_extract_antismash_data(
    url: str, antismash_id: str, download_root: str | PathLike, extract_root: str | PathLike
) -> None:
    """Download and extract antiSMASH BGC archive for a specified genome.

    This function downloads a BGC archive from the specified URL, extracts its contents,
    and organizes the extracted files into a structured directory under the given `extract_root`.

    Args:
        url (str): The URL to download the BGC archive from.
        antismash_id (str): The identifier for the antiSMASH genome, used to name the extraction directory.
        download_root: Path to the directory where the downloaded archive will be stored.
        extract_root: Path to the directory where the data files will be extracted.
            Note that an `antismash` directory will be created in the specified `extract_root` if
            it doesn't exist. The files will be extracted to `<extract_root>/antismash/<antismash_id>` directory.

    Raises:
        ValueError: if `<extract_root>/antismash/<antismash_id>` dir is not empty.
        Exception: If any error occurs during the download or extraction process, the partially extracted
            directory will be cleaned up, and the exception will be re-raised.

    Examples:
         >>> download_and_extract_antismash_data(
                 "https://antismash-db.secondarymetabolites.org/output/GCF_001.1/GCF_001.1.zip",
                 "GCF_001.1",
                 "/data/download",
                 "/data/extracted"
             )
    """
    extract_path = Path(extract_root) / "antismash" / antismash_id

    _prepare_extract_path(extract_path)
    try:
        download_and_extract_archive(url, download_root, extract_path, f"{antismash_id}.zip")
        _cleanup_extracted_files(extract_path)
    except Exception as e:
        shutil.rmtree(extract_path)
        raise e


def download_and_extract_from_antismash_api(
    job_id: str, antismash_id: str, download_root: str | PathLike, extract_root: str | PathLike
) -> None:
    """Downloads and extracts results from an antiSMASH API job.

    This function constructs the download URL using the provided job ID then
    downloads the results as a ZIP file and extracts its contents to the specified directories.

    Args:
        antismash_id (str): The unique identifier for the antiSMASH dataset.
        job_id (str): The job ID for the antiSMASH API job.
        download_root (str or PathLike): The root directory where the ZIP file will be downloaded.
        extract_root (str or PathLike): The root directory where the contents of the ZIP file will be extracted.

    Raises:
        requests.exceptions.RequestException: If there is an issue with the HTTP request.
        zipfile.BadZipFile: If the downloaded file is not a valid ZIP file.
        OSError: If there is an issue with file operations such as writing or extracting.
    """
    url = ANTISMASH_API_DOWNLOAD_URL.format(job_id, antismash_id + ".zip")
    download_and_extract_antismash_data(url, antismash_id, download_root, extract_root)


def download_and_extract_from_antismash_db(
    refseq_acc: str, download_root: str | PathLike, extract_root: str | PathLike
) -> None:
    """Download and extract antiSMASH BGC archive for a specified genome.

    The antiSMASH database (https://antismash-db.secondarymetabolites.org/)
    is used to download the BGC archive. And antiSMASH use RefSeq assembly id
    of a genome as the id of the archive.

    Args:
        refseq_acc: The id used to download BGC archive from antiSMASH database.
            If the id is versioned (e.g., "GCF_004339725.1") please be sure to
            specify the version as well.
        download_root: Path to the directory to place downloaded archive in.
        extract_root: Path to the directory data files will be extracted to.
            Note that an `antismash` directory will be created in the specified `extract_root` if
            it doesn't exist. The files will be extracted to `<extract_root>/antismash/<antismash_id>` directory.

    Raises:
        ValueError: if `<extract_root>/antismash/<refseq_acc>` dir is not empty.

    Examples:
        >>> download_and_extract_from_antismash_db("GCF_004339725.1", "/data/download", "/data/extracted")
    """
    for base_url in [ANTISMASH_DB_DOWNLOAD_URL, ANTISMASH_DBV2_DOWNLOAD_URL]:
        url = base_url.format(refseq_acc, f"{refseq_acc}.zip")
        if requests.head(url).status_code == 404:  # not found
            continue
        download_and_extract_antismash_data(url, refseq_acc, download_root, extract_root)
        return  # Exit the loop once a valid URL is processed

    # if both urls give 404 not found
    raise RuntimeError(f"No results in antiSMASH DB for {refseq_acc}")


def _check_extract_path(extract_path: Path):
    # check if extract_path is empty
    if any(extract_path.iterdir()):
        raise ValueError(f'Nonempty directory: "{extract_path}"')


def _cleanup_extracted_files(extract_path: str | PathLike) -> None:
    # delete subdirs
    for subdir_path in list_dirs(extract_path):
        shutil.rmtree(subdir_path)

    # delete unnecessary files
    files_to_keep = list_files(extract_path, suffix=(".json", ".gbk"))
    for file in list_files(extract_path):
        if file not in files_to_keep:
            os.remove(file)


def _prepare_extract_path(extract_path: str | PathLike) -> None:
    if extract_path.exists():
        _check_extract_path(extract_path)
    else:
        extract_path.mkdir(parents=True, exist_ok=True)
