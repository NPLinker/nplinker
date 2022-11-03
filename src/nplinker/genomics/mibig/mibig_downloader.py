import os
import shutil
from nplinker.logconfig import LogConfig
from nplinker.utils import download_and_extract_archive

logger = LogConfig.getLogger(__file__)

MIBIG_METADATA_URL = "https://dl.secondarymetabolites.org/mibig/mibig_json_{version}.tar.gz"

_MD5_MIBIG_METADATA = {
    "1.0": "035a14e94d2733eb61f615f418c08494",
    "1.1": "e63ceca82363ac27d50650e133ae3fa1",
    "1.2": "46e862018bd076d0b6072e62b7d8cfa2",
    "1.3": "94b3a761323709d06b663232f30210b4",
    "1.4": "a85530571d9dd7978b1bb0f2580cd30e",
    "2.0": "843ce4677db6d11422f0e6d94dd03e81",
    "3.0": "7c38b90f939086c03392d99a913baef9",
    "3.1": "643d1349722a9437d8dcf558dac5f815"
}


def download_and_extract_mibig_metadata(
    download_root: str,
    extract_path: str,
    version: str = "3.1",
):
    """Download and extract MIBiG metadata json files

    Args:
        download_root(str): Path to the directory to place downloaded archive in
        extract_path(str): Path to the empty directory metadata json files will
            be extracted to.
            Note that if it's a directory containing the metadata json files
            after extraction, the json files will be moved to `extract_path`,
            and then the extracted directory will be removed.
        version (str, optional): _description_. Defaults to "3.1".

    Examples:
        >>> download_and_extract_mibig_metadata("/data/download", "/data/mibig_metadata")
    """
    if download_root == extract_path:
        raise ValueError(
            "Identical path of download directory and extract directory")

    # check if extract_path is empty (not check dot files)
    files = [i for i in os.listdir(extract_path) if not i.startswith(".")]
    if len(files) != 0:
        raise ValueError(f'Nonempty directory: "{extract_path}"')

    # download and extract
    md5 = _MD5_MIBIG_METADATA[version]
    download_and_extract_archive(
        url=MIBIG_METADATA_URL.format(version=version),
        download_root=download_root,
        extract_root=extract_path,
        md5=md5)

    # After extracting mibig archive, it's either one dir or many json files,
    # if it's a dir, then move all json files from it to extract_path
    subdirs = [i for i in os.listdir(extract_path) if not i.startswith(".")]
    if len(subdirs) == 1:
        subdir_name = subdirs[0]
        subdir_path = os.path.join(extract_path, subdir_name)
        for fname in os.listdir(subdir_path):
            if fname.startswith("BGC"):  # to avoid dot files
                shutil.move(os.path.join(subdir_path, fname),
                            os.path.join(extract_path, fname))

        # delete subdir
        if subdir_path != extract_path:
            shutil.rmtree(subdir_path)
