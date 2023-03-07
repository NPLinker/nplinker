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

# TODO: add unit tests (only public func)
# TODO: add doc string
def download_and_extract_antismash_metadata(
        refseq_assembly_id: str,
        download_root: str,
        extract_path: str):
    """_summary_

    Args:
        refseq_assembly_id(str): _description_
        download_root(str): _description_
        extract_path(str): _description_

    Raises:
        ValueError: _description_
        ValueError: _description_

    Examples:
        >>> 
        """

    if not os.path.exists(extract_path):
        os.makedirs(extract_path, exist_ok=True)

    if download_root == extract_path:
        raise ValueError(
            "Identical path of download directory and extract directory")

    # check if extract_path is empty
    files = [i for i in os.listdir(extract_path)]
    if len(files) != 0:
        raise ValueError(f'Nonempty directory: "{extract_path}"')

    for base_url in [
            ANTISMASH_DB_DOWNLOAD_URL, ANTISMASH_DBV2_DOWNLOAD_URL
    ]:
        url = base_url.format(refseq_assembly_id, refseq_assembly_id + '.zip')

        download_and_extract_archive(url, download_root, extract_path, refseq_assembly_id + '.zip')
        logger.info(
            f'Genome data successfully extracted for {refseq_assembly_id}')
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
    logger.info(f'download_and_extract_antismash_metadata process for {refseq_assembly_id} is over')
