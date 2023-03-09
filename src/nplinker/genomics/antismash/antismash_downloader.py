import os
import shutil
from nplinker.logconfig import LogConfig
from nplinker.utils import download_and_extract_archive, list_dirs, list_files

logger = LogConfig.getLogger(__name__)

# urls to be given to download antismash data
ANTISMASH_DB_PAGE_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/'
ANTISMASH_DB_DOWNLOAD_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/{}'

# The antiSMASH DBV2 is for the availability of the old version, better to keep it.
ANTISMASH_DBV2_PAGE_URL = 'https://antismash-dbv2.secondarymetabolites.org/output/{}/'
ANTISMASH_DBV2_DOWNLOAD_URL = 'https://antismash-dbv2.secondarymetabolites.org/output/{}/{}'


def _check_roots(download_root, extract_root):
    if download_root == extract_root:
        raise ValueError(
            "Identical path of download directory and extract directory")


def _check_extract_path(extract_path):
    if os.path.exists(extract_path):
        # check if extract_path is empty
        files = list(os.listdir(extract_path))
        if len(files) != 0:
            raise ValueError(f'Nonempty directory: "{extract_path}"')
    else:
        os.makedirs(extract_path, exist_ok=True)


def download_and_extract_antismash_metadata(refseq_assembly_id: str,
                                            download_root: str,
                                            extract_root: str):
    """Download and extract Antismash files for a specified refseq_assembly_id.

    Args:
        refseq_assembly_id(str): Assembled genomic RefSeq (reference sequence) id.
            If the id is versioned (e.g., "GCF_004339725.1") please be sure to 
            specify the version as well. 
        download_root(str): Path to the directory to place downloaded archive in.
        extract_root(str): Path to the directory data files will be extracted to.
            Note that if it will create an antismash/ directory in the specified extract_root, if
            it doesn't already exist.
            The files will be extracted to <extract_root>/antismash/<refseq_assembly_id>/ dir.

    Raises:
        ValueError: if download_root and extract_root dirs are the same.
        ValueError: if <extract_root>/antismash/<refseq_assembly_id> dir is not empty.

    Examples:
        >>> download_and_extract_antismash_metadata("GCF_004339725.1", "/data/download", "/data/extracted")
        """
    extract_path = os.path.join(extract_root, "antismash", refseq_assembly_id)

    _check_roots(download_root, extract_root)
    _check_extract_path(extract_path)

    for base_url in [ANTISMASH_DB_DOWNLOAD_URL, ANTISMASH_DBV2_DOWNLOAD_URL]:
        url = base_url.format(refseq_assembly_id, refseq_assembly_id + '.zip')

        download_and_extract_archive(url, download_root, extract_path,
                                     refseq_assembly_id + '.zip')
        logger.info(
            'Genome data successfully extracted for %s', refseq_assembly_id)
        break

    # delete subdirs
    logger.info('Deleting unnecessary extracted subdirs and files')
    subdirs = list_dirs(extract_path)
    for subdir_path in subdirs:
        shutil.rmtree(subdir_path)

    files_to_keep = list_files(extract_path, suffix=(".json", ".gbk"))

    for file in list_files(extract_path):
        if file not in files_to_keep:
            os.remove(file)
    logger.info(
        'download_and_extract_antismash_metadata process for %s is over', refseq_assembly_id
    )
