import ast
import csv
import os
import re
import time
import zipfile
from os import PathLike
from typing import Dict
import httpx
from bs4 import BeautifulSoup
from bs4 import NavigableString
from bs4 import Tag
from deprecated import deprecated
from progress.bar import Bar
from nplinker.logconfig import LogConfig
from . import download_and_extract_antismash_data

logger = LogConfig.getLogger(__name__)

# urls to be given to download antismash data
ANTISMASH_DB_PAGE_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/'
ANTISMASH_DB_DOWNLOAD_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/{}'

# The antiSMASH DBV2 is for the availability of the old version, better to keep it.
ANTISMASH_DBV2_PAGE_URL = 'https://antismash-dbv2.secondarymetabolites.org/output/{}/'
ANTISMASH_DBV2_DOWNLOAD_URL = 'https://antismash-dbv2.secondarymetabolites.org/output/{}/{}'

NCBI_LOOKUP_URL_NEW = 'https://www.ncbi.nlm.nih.gov/assembly/?term={}'

JGI_GENOME_LOOKUP_URL = 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid={}'

USER_AGENT = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0'


class GenomeStatus:

    def __init__(self,
                 original_id,
                 resolved_refseq_id,
                 attempted=False,
                 filename=""):
        self.original_id = ';'.join(original_id.split(','))
        self.resolved_refseq_id = None if resolved_refseq_id == 'None' else resolved_refseq_id
        if ast.literal_eval(attempted):
            self.attempted = True
        else:
            self.attempted = False
        self.filename = filename

    @classmethod
    def from_csv(cls, original_id, resolved_refseq_id, attempted, filename):
        return cls(original_id, resolved_refseq_id, attempted, filename)

    def to_csv(self):
        return ','.join([
            str(self.original_id),
            str(self.resolved_refseq_id),
            str(self.attempted), self.filename
        ])


def _get_genome_status_log(
        genome_status_file: str | PathLike) -> Dict[str, GenomeStatus]:
    """Read genome_status.txt in a dict if it exists, otherwise create
    the dictionary. 

    Args:
        genome_status_file(str | PathLike): it records genome IDs and local filenames to avoid having to repeat HTTP requests
        each time the app is loaded (this can take a lot of time if there are dozens of genomes).

    Returns:
        dict: log for information about the read genome IDs.
        """

    genome_status = {}

    # 'genome_status.txt' is read, then in the for loop over the genome records it gets updated,
    # and finally it is saved again in 'genome_status.txt' which is overwritten
    if os.path.exists(genome_status_file):
        with open(genome_status_file) as f:
            for line in csv.reader(f):
                asobj = GenomeStatus.from_csv(*line)
                genome_status[asobj.original_id] = asobj

    return genome_status


def _get_best_available_genome_id(
        genome_id_data: Dict[str, str]) -> str | None:
    """Get the best available ID from genome_id_data dict.

    Args:
        genome_id_data(dict): dictionary containing information
        for each genome record present.

    Returns:
        str | None: ID for the genome, if present, otherwise None. 

        """
    if 'RefSeq_accession' in genome_id_data:
        return genome_id_data['RefSeq_accession']
    if 'GenBank_accession' in genome_id_data:
        return genome_id_data['GenBank_accession']
    if 'JGI_Genome_ID' in genome_id_data:
        return genome_id_data['JGI_Genome_ID']

    logger.warning(
        f'No known genome ID field in genome data: {genome_id_data}')
    return None


def _ncbi_genbank_search(
        genbank_id: str,
        retry_time: float = 5.0) -> Tag | NavigableString | None:
    url = NCBI_LOOKUP_URL_NEW.format(genbank_id)
    logger.debug(f'Looking up GenBank data for {genbank_id} at {url}')
    resp = httpx.get(url, follow_redirects=True)

    if resp.status_code == httpx.codes.OK:
        # the page should contain a <dl> element with class "assembly_summary_new". retrieving
        # the page seems to fail occasionally in the middle of lengthy sequences of genome
        # lookups, so there might be some throttling going on. this will automatically retry
        # the lookup if the expected content isn't found the first time
        soup = BeautifulSoup(resp.content, 'html.parser')
        # find the <dl> element with class "assembly_summary_new"
        dl_element = soup.find('dl', {'class': 'assembly_summary_new'})
        if dl_element is not None:
            return dl_element

    logger.debug(
        f'NCBI lookup failed, status code {resp.status_code}. Trying again in {retry_time} seconds'
    )
    time.sleep(retry_time)
    logger.debug(f'Looking up GenBank data for {genbank_id} at {url}')
    resp = httpx.get(url, follow_redirects=True)
    if resp.status_code == httpx.codes.OK:
        soup = BeautifulSoup(resp.content, 'html.parser')
        # find the <dl> element with class "assembly_summary_new"
        dl_element = soup.find('dl', {'class': 'assembly_summary_new'})
        if dl_element is not None:
            return dl_element

    logger.warning(
        f'Failed to resolve NCBI genome ID {genbank_id} at URL {url} (after retrying)'
    )
    return None


def _resolve_genbank_accession(genbank_id: str) -> str | None:
    """Try to get RefSeq id through given GenBank id. 

    Args:
        genbank_id(str): ID for GenBank accession. 

    Raises:
        Exception: "Unknown HTML format" if the search of genbank does not give any results. 
        Exception: "Expected HTML elements not found" if no match with RefSeq assembly accession is found.

    Returns:
        str | None: RefSeq ID if the search is successful, otherwise None. 
        """
    logger.info(
        f'Attempting to resolve Genbank accession {genbank_id} to RefSeq accession'
    )
    # genbank id => genbank seq => refseq

    # The GenBank accession can have several formats:
    # 1: BAFR00000000.1
    # 2: NZ_BAGG00000000.1
    # 3: NC_016887.1
    # Case 1 is the default.
    if '_' in genbank_id:
        # case 2
        if len(genbank_id.split('_')[-1].split('.')[0]) == 12:
            genbank_id = genbank_id.split('_')[-1]
        # case 3
        else:
            genbank_id = genbank_id.lower()

    # get rid of any extraneous whitespace
    genbank_id = genbank_id.strip()
    logger.debug(f'Parsed GenBank ID to "{genbank_id}"')

    # run a search using the GenBank accession ID
    try:
        dl_element = _ncbi_genbank_search(genbank_id)
        if dl_element is None:
            raise Exception('Unknown HTML format')

        refseq_idx = -1
        for field_idx, field in enumerate(dl_element.children):
            # this is the element immediately preceding the one with
            # the actual RefSeq ID we want
            if field.getText().strip() == 'RefSeq assembly accession:':
                refseq_idx = field_idx + 1

            # this should be True when we've reached the right element
            if field_idx == refseq_idx:
                refseq_id = field.getText()
                # if it has any spaces, take everything up to first one (some have annotations afterwards)
                if refseq_id.find(' ') != -1:
                    refseq_id = refseq_id[:refseq_id.find(' ')]

                return refseq_id

        if refseq_idx == -1:
            raise Exception('Expected HTML elements not found')
    except Exception as e:
        logger.warning(
            f'Failed resolving GenBank accession {genbank_id}, error {e}')

    return None


def _resolve_jgi_accession(jgi_id: str) -> str | None:
    """Try to get RefSeq id through given JGI id.

    Args:
        jgi_id(str): JGI_Genome_ID for GenBank accession. 

    Returns:
        str | None: GenBank accession ID if search is successful, otherwise None.
        """
    url = JGI_GENOME_LOOKUP_URL.format(jgi_id)
    logger.info(
        f'Attempting to resolve JGI_Genome_ID {jgi_id} to GenBank accession via {url}'
    )
    # no User-Agent header produces a 403 Forbidden error on this site...
    try:
        resp = httpx.get(url,
                         headers={'User-Agent': USER_AGENT},
                         timeout=10.0,
                         follow_redirects=True)
    except httpx.ReadTimeout:
        logger.warning('Timed out waiting for result of JGI_Genome_ID lookup')
        return None

    soup = BeautifulSoup(resp.content, 'html.parser')
    # find the table entry giving the NCBI assembly accession ID
    link = soup.find(
        'a', href=re.compile('https://www.ncbi.nlm.nih.gov/nuccore/.*'))
    if link is None:
        return None

    return _resolve_genbank_accession(link.text)


def _resolve_refseq_id(genome_id_data: Dict[str, str]) -> str | None:
    """Get the RefSeq ID to which the genome accession is linked.
    Check https://pairedomicsdata.bioinformatics.nl/schema.json.

    Args:
        genome_id_data(dict): dictionary containing information
        for each genome record present.

    Returns:
        str | None: ID for the accession genome, if present, otherwise None. 
        """
    if 'RefSeq_accession' in genome_id_data:
        # best case, can use this directly
        return genome_id_data['RefSeq_accession']
    if 'GenBank_accession' in genome_id_data:
        # resolve via NCBI
        return _resolve_genbank_accession(genome_id_data['GenBank_accession'])
    if 'JGI_Genome_ID' in genome_id_data:
        # resolve via JGI => NCBI
        return _resolve_jgi_accession(genome_id_data['JGI_Genome_ID'])

    logger.warning(f'Unable to resolve genome_ID: {genome_id_data}')
    return None


def _get_antismash_db_page(genome_obj: GenomeStatus) -> str | None:
    # want to try up to 4 different links here, v1 and v2 databases, each
    # with and without the .1 suffix on the accesssion ID

    accesssions = [
        genome_obj.resolved_refseq_id, genome_obj.resolved_refseq_id + '.1'
    ]
    for base_url in [ANTISMASH_DB_PAGE_URL, ANTISMASH_DBV2_PAGE_URL]:
        for accession in accesssions:
            url = base_url.format(accession)
            link = None

            logger.info(f'antismash DB lookup for {accession} => {url}')
            try:
                resp = httpx.get(url, follow_redirects=True)
                soup = BeautifulSoup(resp.content, 'html.parser')
                # retrieve .zip file download link
                link = soup.find('a',
                                 {'href': lambda url: url.endswith('.zip')})
            except Exception as e:
                logger.debug(f'antiSMASH DB page load failed: {e}')

            if link is not None:
                logger.info(
                    f"antiSMASH lookup succeeded! Filename is {link['href']}")
                # save with the .1 suffix if that worked
                genome_obj.resolved_refseq_id = accession
                return link['href']

    return None


def _get_antismash_zip_data(accession_id: str, filename: str,
                            local_path: str | PathLike) -> bool:
    for base_url in [ANTISMASH_DB_DOWNLOAD_URL, ANTISMASH_DBV2_DOWNLOAD_URL]:
        zipfile_url = base_url.format(accession_id, filename)
        with open(local_path, 'wb') as f:
            total_bytes = 0
            try:
                with httpx.stream('GET', zipfile_url) as r:
                    if r.status_code == 404:
                        logger.debug('antiSMASH download URL was a 404')
                        continue

                    logger.info(f'Downloading from antiSMASH: {zipfile_url}')
                    filesize = int(r.headers['content-length'])
                    bar_obj = Bar(filename,
                                  max=filesize,
                                  suffix='%(percent)d%%')
                    for data in r.iter_bytes():
                        f.write(data)
                        total_bytes += len(data)
                        bar_obj.next(len(data))
                    bar_obj.finish()
            except Exception as e:
                logger.warning(f'antiSMASH zip download failed: {e}')
                continue

        return True

    return False


def _download_antismash_zip(antismash_obj: GenomeStatus,
                            project_download_cache: str | PathLike) -> bool:
    # save zip files to avoid having to repeat above lookup every time
    local_path = os.path.join(project_download_cache,
                              f'{antismash_obj.resolved_refseq_id}.zip')
    logger.debug(f'Checking for existing antismash zip at {local_path}')

    cached = False
    # if the file exists locally
    if os.path.exists(local_path):
        logger.info(f'Found cached file at {local_path}')
        try:
            # check if it's a valid zip file, if so treat it as cached
            with zipfile.ZipFile(local_path) as _:
                cached = True
                antismash_obj.filename = local_path
        except zipfile.BadZipFile as bzf:
            # otherwise delete and redownload
            logger.info(
                f'Invalid antismash zipfile found ({bzf}). Will download again'
            )
            os.unlink(local_path)
            antismash_obj.filename = ""

    if not cached:
        filename = _get_antismash_db_page(antismash_obj)
        if filename is None:
            return False

        _get_antismash_zip_data(antismash_obj.resolved_refseq_id, filename,
                                local_path)
        antismash_obj.filename = local_path

    return True


def _extract_antismash_zip(antismash_obj: GenomeStatus,
                           project_file_cache: str | PathLike) -> bool:
    if antismash_obj.filename is None or len(antismash_obj.filename) == 0:
        return False

    output_path = os.path.join(project_file_cache, 'antismash',
                               antismash_obj.resolved_refseq_id)
    exists_already = os.path.exists(output_path) and os.path.exists(
        os.path.join(output_path, 'completed'))

    logger.debug(
        f'Extracting antismash data to {output_path}, exists_already = {exists_already}'
    )
    if exists_already:
        return True

    # create a subfolder for each set of genome data (the zip files used to be
    # constructed with path info but that seems to have changed recently)
    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)

    with zipfile.ZipFile(antismash_obj.filename) as antismash_zip:
        kc_prefix1 = f'{antismash_obj.resolved_refseq_id}/knownclusterblast'
        kc_prefix2 = 'knownclusterblast'
        for zip_member in antismash_zip.namelist():
            # TODO other files here?
            if zip_member.endswith('.gbk') or zip_member.endswith('.json'):
                antismash_zip.extract(zip_member, path=output_path)
            elif zip_member.startswith(kc_prefix1) or zip_member.startswith(
                    kc_prefix2):
                if zip_member.endswith(
                        '.txt') and 'mibig_hits' not in zip_member:
                    antismash_zip.extract(zip_member, path=output_path)

    with open(os.path.join(output_path, 'completed'), 'w'):
        pass

    return True


def podp_download_and_extract_antismash_data(
        genome_records: list[Dict[str, Dict[str, str]]],
        project_download_cache: str | PathLike,
        project_file_cache: str | PathLike):

    genome_status_file = os.path.join(project_download_cache,
                                      'genome_status.txt')
    genome_status = _get_genome_status_log(genome_status_file)

    for i, genome_record in enumerate(genome_records):
        # get the best available ID from the dict
        raw_genome_id = _get_best_available_genome_id(
            genome_record['genome_ID'])
        if raw_genome_id is None:
            logger.warning(
                f'Ignoring genome record "{genome_record}" due to missing genome ID field'
            )
            continue

        # use this to check if the lookup has already been attempted and if
        # so if the file is cached locally
        if raw_genome_id not in genome_status:
            genome_status[raw_genome_id] = GenomeStatus(raw_genome_id, None)

        genome_obj = genome_status[raw_genome_id]

        logger.info(
            f'Checking for antismash data {i + 1}/{len(genome_records)}, current genome ID={raw_genome_id}'
        )
        # first check if file is cached locally
        if os.path.exists(genome_obj.filename):
            # file already downloaded
            logger.info(
                f'Genome ID {raw_genome_id} already downloaded to {genome_obj.filename}'
            )
            genome_record['resolved_refseq_id'] = genome_obj.resolved_refseq_id
        elif genome_obj.attempted:
            # lookup attempted previously but failed
            logger.info(
                f'Genome ID {raw_genome_id} skipped due to previous failure')
            genome_record['resolved_refseq_id'] = genome_obj.resolved_refseq_id
        else:
            # if no existing file and no lookup attempted, can start process of
            # trying to retrieve the data

            # lookup the ID
            logger.info(
                f'Beginning lookup process for genome ID {raw_genome_id}')

            genome_obj.resolved_refseq_id = _resolve_refseq_id(
                genome_record['genome_ID'])
            genome_obj.attempted = True

            if genome_obj.resolved_refseq_id is None:
                # give up on this one
                logger.warning(f'Failed lookup for genome ID {raw_genome_id}')
                with open(genome_status_file, 'a+') as f:
                    f.write(genome_obj.to_csv() + '\n')
                continue

        # if we got a refseq ID, now try to download and extract the data from antismash
        download_and_extract_antismash_data(genome_obj.resolved_refseq_id,
                                            project_download_cache,
                                            project_file_cache)

        with open(genome_status_file, 'a+', newline='\n') as f:
            f.write(genome_obj.to_csv() + '\n')

        output_path = os.path.join(project_file_cache, 'antismash',
                                   genome_obj.resolved_refseq_id)
        if not os.path.exists(os.path.join(output_path, 'completed')):
            with open(os.path.join(output_path, 'completed'), 'w'):
                pass

    missing = len([x for x in genome_status.values() if len(x.filename) == 0])
    logger.info(
        f'Dataset has {missing} missing sets of antiSMASH data (from a total of {len(genome_records)})'
    )

    with open(genome_status_file, 'w', newline='\n') as f:
        for obj in genome_status.values():
            f.write(obj.to_csv() + '\n')

    if missing == len(genome_records):
        logger.warning('Failed to successfully retrieve ANY genome data!')


@deprecated(version="1.3.3",
            reason="Use download_and_extract_antismash_data class instead.")
def download_antismash_data(genome_records: list[Dict[str, Dict[str, str]]],
                            project_download_cache: str | PathLike,
                            project_file_cache: str | PathLike):

    genome_status_file = os.path.join(project_download_cache,
                                      'genome_status.txt')
    genome_status = _get_genome_status_log(genome_status_file)

    for i, genome_record in enumerate(genome_records):
        # get the best available ID from the dict
        raw_genome_id = _get_best_available_genome_id(
            genome_record['genome_ID'])
        if raw_genome_id is None:
            logger.warning(
                f'Ignoring genome record "{genome_record}" due to missing genome ID field'
            )
            continue

        # use this to check if the lookup has already been attempted and if
        # so if the file is cached locally
        if raw_genome_id not in genome_status:
            genome_status[raw_genome_id] = GenomeStatus(raw_genome_id, None)

        genome_obj = genome_status[raw_genome_id]

        logger.info(
            f'Checking for antismash data {i + 1}/{len(genome_records)}, current genome ID={raw_genome_id}'
        )
        # first check if file is cached locally
        if os.path.exists(genome_obj.filename):
            # file already downloaded
            logger.info(
                f'Genome ID {raw_genome_id} already downloaded to {genome_obj.filename}'
            )
            genome_record['resolved_refseq_id'] = genome_obj.resolved_refseq_id
        elif genome_obj.attempted:
            # lookup attempted previously but failed
            logger.info(
                f'Genome ID {raw_genome_id} skipped due to previous failure')
            genome_record['resolved_refseq_id'] = genome_obj.resolved_refseq_id
        else:
            # if no existing file and no lookup attempted, can start process of
            # trying to retrieve the data

            # lookup the ID
            logger.info(
                f'Beginning lookup process for genome ID {raw_genome_id}')

            genome_obj.resolved_refseq_id = _resolve_refseq_id(
                genome_record['genome_ID'])
            genome_obj.attempted = True

            if genome_obj.resolved_refseq_id is None:
                # give up on this one
                logger.warning(f'Failed lookup for genome ID {raw_genome_id}')
                with open(genome_status_file, 'a+') as f:
                    f.write(genome_obj.to_csv() + '\n')
                continue

            # if we got a refseq ID, now try to download the data from antismash
            if _download_antismash_zip(genome_obj, project_download_cache):
                logger.info(
                    f'Genome data successfully downloaded for {raw_genome_id}')
                genome_record[
                    'resolved_refseq_id'] = genome_obj.resolved_refseq_id
            else:
                logger.warning(
                    f'Failed to download antiSMASH data for genome ID {genome_obj.resolved_refseq_id} ({genome_obj.original_id})'
                )

            with open(genome_status_file, 'a+', newline='\n') as f:
                f.write(genome_obj.to_csv() + '\n')

        _extract_antismash_zip(genome_obj, project_file_cache)

    missing = len([x for x in genome_status.values() if len(x.filename) == 0])
    logger.info(
        f'Dataset has {missing} missing sets of antiSMASH data (from a total of {len(genome_records)})'
    )

    with open(genome_status_file, 'w', newline='\n') as f:
        for obj in genome_status.values():
            f.write(obj.to_csv() + '\n')

    if missing == len(genome_records):
        logger.warning('Failed to successfully retrieve ANY genome data!')
