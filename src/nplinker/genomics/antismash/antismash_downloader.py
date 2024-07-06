from __future__ import annotations
import logging
import os
import shutil
from os import PathLike
from pathlib import Path
from nplinker.utils import download_and_extract_archive
from nplinker.utils import list_dirs
from nplinker.utils import list_files


logger = logging.getLogger(__name__)

# urls to be given to download antismash data
ANTISMASH_DB_DOWNLOAD_URL = "https://antismash-db.secondarymetabolites.org/output/{}/{}"
# The antiSMASH DBV2 is for the availability of the old version, better to keep it.
ANTISMASH_DBV2_DOWNLOAD_URL = "https://antismash-dbv2.secondarymetabolites.org/output/{}/{}"


def download_and_extract_antismash_data(
    antismash_id: str, download_root: str | PathLike, extract_root: str | PathLike
) -> None:
    """Download and extract antiSMASH BGC archive for a specified genome.

    The antiSMASH database (https://antismash-db.secondarymetabolites.org/)
    is used to download the BGC archive. And antiSMASH use RefSeq assembly id
    of a genome as the id of the archive.

    Args:
        antismash_id: The id used to download BGC archive from antiSMASH database.
            If the id is versioned (e.g., "GCF_004339725.1") please be sure to
            specify the version as well.
        download_root: Path to the directory to place downloaded archive in.
        extract_root: Path to the directory data files will be extracted to.
            Note that an `antismash` directory will be created in the specified `extract_root` if
            it doesn't exist. The files will be extracted to `<extract_root>/antismash/<antismash_id>` directory.

    Raises:
        ValueError: if `<extract_root>/antismash/<refseq_assembly_id>` dir is not empty.

    Examples:
        >>> download_and_extract_antismash_metadata("GCF_004339725.1", "/data/download", "/data/extracted")
    """
    download_root = Path(download_root)
    extract_root = Path(extract_root)
    extract_path = extract_root / "antismash" / antismash_id

    try:
        if extract_path.exists():
            _check_extract_path(extract_path)
        else:
            extract_path.mkdir(parents=True, exist_ok=True)

        for base_url in [ANTISMASH_DB_DOWNLOAD_URL, ANTISMASH_DBV2_DOWNLOAD_URL]:
            url = base_url.format(antismash_id, antismash_id + ".zip")
            download_and_extract_archive(url, download_root, extract_path, antismash_id + ".zip")
            break

        # delete subdirs
        for subdir_path in list_dirs(extract_path):
            shutil.rmtree(subdir_path)

        # delete unnecessary files
        files_to_keep = list_files(extract_path, suffix=(".json", ".gbk"))
        for file in list_files(extract_path):
            if file not in files_to_keep:
                os.remove(file)

        logger.info("antiSMASH BGC data of %s is downloaded and extracted.", antismash_id)

    except Exception as e:
        shutil.rmtree(extract_path)
        logger.warning(e)
        raise e


def _check_extract_path(extract_path: Path):
    # check if extract_path is empty
    if any(extract_path.iterdir()):
        raise ValueError(f'Nonempty directory: "{extract_path}"')
