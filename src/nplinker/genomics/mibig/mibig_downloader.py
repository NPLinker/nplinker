from __future__ import annotations
import logging
import os
import shutil
from pathlib import Path
from nplinker.utils import download_and_extract_archive
from nplinker.utils import list_dirs
from nplinker.utils import list_files


logger = logging.getLogger(__name__)

MIBIG_METADATA_URL = "https://dl.secondarymetabolites.org/mibig/mibig_json_{version}.tar.gz"

_MD5_MIBIG_METADATA = {
    "1.0": "035a14e94d2733eb61f615f418c08494",
    "1.1": "e63ceca82363ac27d50650e133ae3fa1",
    "1.2": "46e862018bd076d0b6072e62b7d8cfa2",
    "1.3": "94b3a761323709d06b663232f30210b4",
    "1.4": "a85530571d9dd7978b1bb0f2580cd30e",
    "2.0": "843ce4677db6d11422f0e6d94dd03e81",
    "3.0": "7c38b90f939086c03392d99a913baef9",
    "3.1": "643d1349722a9437d8dcf558dac5f815",
}


def download_and_extract_mibig_metadata(
    download_root: str | os.PathLike,
    extract_path: str | os.PathLike,
    version: str = "3.1",
):
    """Download and extract MIBiG metadata json files.

    Note that it does not matter whether the metadata json files are in nested folders or not in the archive,
    all json files will be extracted to the same location, i.e. `extract_path`. The nested
    folders will be removed if they exist. So the `extract_path` will have only json files.

    Args:
        download_root: Path to the directory in which to place the downloaded archive.
        extract_path: Path to an empty directory where the json files will be extracted.
            The directory must be empty if it exists. If it doesn't exist, the directory will be created.
        version: _description_. Defaults to "3.1".

    Examples:
        >>> download_and_extract_mibig_metadata("/data/download", "/data/mibig_metadata")
    """
    download_root = Path(download_root)
    extract_path = Path(extract_path)

    if download_root == extract_path:
        raise ValueError("Identical path of download directory and extract directory")

    # check if extract_path is empty
    if not extract_path.exists():
        extract_path.mkdir(parents=True)
    else:
        if len(list(extract_path.iterdir())) != 0:
            raise ValueError(f'Nonempty directory: "{extract_path}"')

    # download and extract
    md5 = _MD5_MIBIG_METADATA[version]
    download_and_extract_archive(
        url=MIBIG_METADATA_URL.format(version=version),
        download_root=download_root,
        extract_root=extract_path,
        md5=md5,
    )

    # After extracting mibig archive, it's either one dir or many json files,
    # if it's a dir, then move all json files from it to extract_path
    subdirs = list_dirs(extract_path)
    if len(subdirs) > 1:
        raise ValueError(f"Expected one extracted directory, got {len(subdirs)}")

    if len(subdirs) == 1:
        subdir_path = subdirs[0]
        for fname in list_files(subdir_path, prefix="BGC", suffix=".json", keep_parent=False):
            shutil.move(os.path.join(subdir_path, fname), os.path.join(extract_path, fname))
        # delete subdir
        if subdir_path != extract_path:
            shutil.rmtree(subdir_path)
