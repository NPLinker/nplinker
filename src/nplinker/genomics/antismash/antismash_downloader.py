import csv
import os
import re
import time
import zipfile
import httpx
from deprecated import deprecated
from bs4 import BeautifulSoup
from progress.bar import Bar
from nplinker.logconfig import LogConfig


logger = LogConfig.getLogger(__name__)

# urls to be given to download antismash data
ANTISMASH_DB_PAGE_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/'
ANTISMASH_DB_DOWNLOAD_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/{}'

ANTISMASH_DBV2_PAGE_URL = 'https://antismash-dbv2.secondarymetabolites.org/output/{}/'
ANTISMASH_DBV2_DOWNLOAD_URL = 'https://antismash-dbv2.secondarymetabolites.org/output/{}/{}'

NCBI_LOOKUP_URL_NEW = 'https://www.ncbi.nlm.nih.gov/assembly/?term={}'

JGI_GENOME_LOOKUP_URL = 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid={}'

USER_AGENT = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0'

class GenomeStatus:

    def __init__(self, original_id, resolved_id, attempted=False, filename=""):
        self.original_id = ';'.join(original_id.split(','))
        self.resolved_id = None if resolved_id == 'None' else resolved_id
        self.attempted = True if attempted == 'True' else False
        self.filename = filename

    @classmethod
    def from_csv(cls, original_id, resolved_id, attempted, filename):
        return cls(original_id, resolved_id, attempted, filename)

    def to_csv(self):
        return ','.join([
            str(self.original_id),
            str(self.resolved_id),
            str(self.attempted), self.filename
        ])

def _get_best_available_genome_id(genome_id_data):
    if 'RefSeq_accession' in genome_id_data:
        return genome_id_data['RefSeq_accession']
    elif 'GenBank_accession' in genome_id_data:
        return genome_id_data['GenBank_accession']
    elif 'JGI_Genome_ID' in genome_id_data:
        return genome_id_data['JGI_Genome_ID']

    logger.warning('No known genome ID field in genome data: {}'.format(
        genome_id_data))
    return None

def _ncbi_genbank_search(genbank_id, retry_time=5.0):
    url = NCBI_LOOKUP_URL_NEW.format(genbank_id)
    logger.debug('Looking up GenBank data for {} at {}'.format(
        genbank_id, url))
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
        'NCBI lookup failed, status code {}. Trying again in {} seconds'.
        format(resp.status_code, retry_time))
    time.sleep(retry_time)
    logger.debug('Looking up GenBank data for {} at {}'.format(
        genbank_id, url))
    resp = httpx.get(url, follow_redirects=True)
    if resp.status_code == httpx.codes.OK:
        soup = BeautifulSoup(resp.content, 'html.parser')
        # find the <dl> element with class "assembly_summary_new"
        dl_element = soup.find('dl', {'class': 'assembly_summary_new'})
        if dl_element is not None:
            return dl_element

    logger.warning(
        'Failed to resolve NCBI genome ID {} at URL {} (after retrying)'.
        format(genbank_id, url))
    return None

def _resolve_genbank_accession(genbank_id):
    logger.info(
        'Attempting to resolve Genbank accession {} to RefSeq accession'.
        format(genbank_id))
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
            'Failed resolving GenBank accession {}, error {}'.format(
                genbank_id, e))

    return None

def _resolve_jgi_accession(jgi_id):
    url = JGI_GENOME_LOOKUP_URL.format(jgi_id)
    logger.info(
        'Attempting to resolve JGI_Genome_ID {} to GenBank accession via {}'
        .format(jgi_id, url))
    # no User-Agent header produces a 403 Forbidden error on this site...
    try:
        resp = httpx.get(url,
                            headers={'User-Agent': USER_AGENT},
                            timeout=10.0,
                            follow_redirects=True)
    except httpx.ReadTimeout:
        logger.warning(
            'Timed out waiting for result of JGI_Genome_ID lookup')
        return None

    soup = BeautifulSoup(resp.content, 'html.parser')
    # find the table entry giving the NCBI assembly accession ID
    link = soup.find(
        'a', href=re.compile('https://www.ncbi.nlm.nih.gov/nuccore/.*'))
    if link is None:
        return None

    return _resolve_genbank_accession(link.text)

def _resolve_genome_id_data(genome_id_data):
    if 'RefSeq_accession' in genome_id_data:
        # best case, can use this directly
        return genome_id_data['RefSeq_accession']
    elif 'GenBank_accession' in genome_id_data:
        # resolve via NCBI
        return _resolve_genbank_accession(
            genome_id_data['GenBank_accession'])
    elif 'JGI_Genome_ID' in genome_id_data:
        # resolve via JGI => NCBI
        return _resolve_jgi_accession(genome_id_data['JGI_Genome_ID'])

    logger.warning(
        f'Unable to resolve genome_ID: {genome_id_data}')
    return None

def _get_antismash_db_page(genome_obj):
    # want to try up to 4 different links here, v1 and v2 databases, each
    # with and without the .1 suffix on the accesssion ID

    accesssions = [genome_obj.resolved_id, genome_obj.resolved_id + '.1']
    for base_url in [ANTISMASH_DB_PAGE_URL, ANTISMASH_DBV2_PAGE_URL]:
        for accession in accesssions:
            url = base_url.format(accession)
            link = None

            logger.info('antismash DB lookup for {} => {}'.format(
                accession, url))
            try:
                resp = httpx.get(url, follow_redirects=True)
                soup = BeautifulSoup(resp.content, 'html.parser')
                # retrieve .zip file download link
                link = soup.find(
                    'a', {'href': lambda url: url.endswith('.zip')})
            except Exception as e:
                logger.debug(f'antiSMASH DB page load failed: {e}')

            if link is not None:
                logger.info(
                    'antiSMASH lookup succeeded! Filename is {}'.format(
                        link['href']))
                # save with the .1 suffix if that worked
                genome_obj.resolved_id = accession
                return link['href']

    return None

def _get_antismash_zip_data(accession_id, filename, local_path):
    for base_url in [
            ANTISMASH_DB_DOWNLOAD_URL, ANTISMASH_DBV2_DOWNLOAD_URL
    ]:
        zipfile_url = base_url.format(accession_id, filename)
        with open(local_path, 'wb') as f:
            total_bytes = 0
            try:
                with httpx.stream('GET', zipfile_url) as r:
                    if r.status_code == 404:
                        logger.debug('antiSMASH download URL was a 404')
                        continue

                    logger.info('Downloading from antiSMASH: {}'.format(
                        zipfile_url))
                    filesize = int(r.headers['content-length'])
                    bar = Bar(filename,
                                max=filesize,
                                suffix='%(percent)d%%')
                    for data in r.iter_bytes():
                        f.write(data)
                        total_bytes += len(data)
                        bar.next(len(data))
                    bar.finish()
            except Exception as e:
                logger.warning(
                    f'antiSMASH zip download failed: {e}')
                continue

        return True

    return False

def _download_antismash_zip(antismash_obj, project_download_cache):
    # save zip files to avoid having to repeat above lookup every time
    local_path = os.path.join(project_download_cache,
                                f'{antismash_obj.resolved_id}.zip')
    logger.debug(
        f'Checking for existing antismash zip at {local_path}')

    cached = False
    # if the file exists locally
    if os.path.exists(local_path):
        logger.info(f'Found cached file at {local_path}')
        try:
            # check if it's a valid zip file, if so treat it as cached
            _ = zipfile.ZipFile(local_path)
            cached = True
            antismash_obj.filename = local_path
        except zipfile.BadZipFile as bzf:
            # otherwise delete and redownload
            logger.info(
                'Invalid antismash zipfile found ({}). Will download again'
                .format(bzf))
            os.unlink(local_path)
            antismash_obj.filename = ""

    if not cached:
        filename = _get_antismash_db_page(antismash_obj)
        if filename is None:
            return False

        _get_antismash_zip_data(antismash_obj.resolved_id, filename,
                                        local_path)
        antismash_obj.filename = local_path

    return True

def _extract_antismash_zip(antismash_obj, project_file_cache):
    if antismash_obj.filename is None or len(antismash_obj.filename) == 0:
        return False

    output_path = os.path.join(project_file_cache, 'antismash',
                                antismash_obj.resolved_id)
    exists_already = os.path.exists(output_path) and os.path.exists(
        os.path.join(output_path, 'completed'))

    logger.debug(
        'Extracting antismash data to {}, exists_already = {}'.format(
            output_path, exists_already))
    if exists_already:
        return True

    # create a subfolder for each set of genome data (the zip files used to be
    # constructed with path info but that seems to have changed recently)
    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)

    antismash_zip = zipfile.ZipFile(antismash_obj.filename)
    kc_prefix1 = f'{antismash_obj.resolved_id}/knownclusterblast'
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

    open(os.path.join(output_path, 'completed'), 'w').close()

    return True

def download_and_extract_antismash_data(item_id, download_root, extract_path): 
    genome_status = {}

    # this file records genome IDs and local filenames to avoid having to repeat HTTP requests
    # each time the app is loaded (this can take a lot of time if there are dozens of genomes)
    genome_status_file = os.path.join(download_root,
                                        'genome_status.txt')

    # genome lookup status info
    if os.path.exists(genome_status_file):
        with open(genome_status_file) as f:
            for line in csv.reader(f):
                asobj = GenomeStatus.from_csv(*line)
                genome_status[asobj.original_id] = asobj

    # use this to check if the lookup has already been attempted and if
    # so if the file is cached locally
    if item_id not in genome_status:
        genome_status[item_id] = GenomeStatus(item_id, None)

    genome_obj = genome_status[item_id]

    logger.info(
        'Checking for antismash data for genome ID={}'.
        format(item_id))
    # first check if file is cached locally
    if os.path.exists(genome_obj.filename):
        # file already downloaded
        logger.info('Genome ID {} already downloaded to {}'.format(
            item_id, genome_obj.filename))
    elif genome_obj.attempted:
        # lookup attempted previously but failed
        logger.info(
            'Genome ID {} skipped due to previous failure'.format(
                item_id))
    else:
        # if no existing file and no lookup attempted, can start process of
        # trying to retrieve the data

        # lookup the ID
        logger.info('Beginning lookup process for genome ID {}'.format(
            item_id))

        genome_obj.resolved_id = item_id # TO CHECK (Cunliang) not sure if this is what we want; in a general case,
        # I don't think we have different possible ids (as in podp json file, for genome_ID nested dicts),
        # so maybe it makes sense to put genome_obj.resolved_id equal to the item_id and only in podp case do the check
        # (through _resolve_genome_id_data, was done here before) outside this function.
        # If this is true, then I think we need to modify GenomeStatus class attributes logic for original_id and resolved_id,
        # which in this way would be the same thing. Then we should modify also the code below, which assumes original_id
        # to be eventually different
        genome_obj.attempted = True

        if genome_obj.resolved_id is None:
            # give up on this one
            logger.warning(
                f'Failed lookup for genome ID {item_id}')
            with open(genome_status_file, 'a+') as f:
                f.write(genome_obj.to_csv() + '\n')

        # if we got a refseq ID, now try to download the data from antismash
        if _download_antismash_zip(genome_obj, download_root):
            logger.info(
                'Genome data successfully downloaded for {}'.format(
                    item_id))
        else:
            logger.warning(
                'Failed to download antiSMASH data for genome ID {} ({})'
                .format(genome_obj.resolved_id,
                        genome_obj.original_id))

        with open(genome_status_file, 'a+', newline='\n') as f:
            f.write(genome_obj.to_csv() + '\n')

    _extract_antismash_zip(genome_obj, extract_path)

    with open(genome_status_file, 'w', newline='\n') as f:
        for obj in genome_status.values():
            f.write(obj.to_csv() + '\n')

@deprecated(version="1.3.3", reason="Use download_and_extract_antismash_data class instead.")
def download_antismash_data(genome_records, project_download_cache, project_file_cache): 
    genome_status = {}

    # this file records genome IDs and local filenames to avoid having to repeat HTTP requests
    # each time the app is loaded (this can take a lot of time if there are dozens of genomes)
    genome_status_file = os.path.join(project_download_cache,
                                        'genome_status.txt')

    # genome lookup status info
    if os.path.exists(genome_status_file):
        with open(genome_status_file) as f:
            for line in csv.reader(f):
                asobj = GenomeStatus.from_csv(*line)
                genome_status[asobj.original_id] = asobj

    for i, genome_record in enumerate(genome_records):
        # get the best available ID from the dict
        best_id = _get_best_available_genome_id(
            genome_record['genome_ID'])
        if best_id is None:
            logger.warning(
                'Ignoring genome record "{}" due to missing genome ID field'
                .format(genome_record))
            continue

        # use this to check if the lookup has already been attempted and if
        # so if the file is cached locally
        if best_id not in genome_status:
            genome_status[best_id] = GenomeStatus(best_id, None)

        genome_obj = genome_status[best_id]

        logger.info(
            'Checking for antismash data {}/{}, current genome ID={}'.
            format(i + 1, len(genome_records), best_id))
        # first check if file is cached locally
        if os.path.exists(genome_obj.filename):
            # file already downloaded
            logger.info('Genome ID {} already downloaded to {}'.format(
                best_id, genome_obj.filename))
            genome_record['resolved_id'] = genome_obj.resolved_id
        elif genome_obj.attempted:
            # lookup attempted previously but failed
            logger.info(
                'Genome ID {} skipped due to previous failure'.format(
                    best_id))
            genome_record['resolved_id'] = genome_obj.resolved_id
        else:
            # if no existing file and no lookup attempted, can start process of
            # trying to retrieve the data

            # lookup the ID
            logger.info('Beginning lookup process for genome ID {}'.format(
                best_id))

            genome_obj.resolved_id = _resolve_genome_id_data(
                genome_record['genome_ID'])
            genome_obj.attempted = True

            if genome_obj.resolved_id is None:
                # give up on this one
                logger.warning(
                    f'Failed lookup for genome ID {best_id}')
                with open(genome_status_file, 'a+') as f:
                    f.write(genome_obj.to_csv() + '\n')
                continue

            # if we got a refseq ID, now try to download the data from antismash
            if _download_antismash_zip(genome_obj, project_download_cache):
                logger.info(
                    'Genome data successfully downloaded for {}'.format(
                        best_id))
                genome_record['resolved_id'] = genome_obj.resolved_id
            else:
                logger.warning(
                    'Failed to download antiSMASH data for genome ID {} ({})'
                    .format(genome_obj.resolved_id,
                            genome_obj.original_id))

            with open(genome_status_file, 'a+', newline='\n') as f:
                f.write(genome_obj.to_csv() + '\n')

        _extract_antismash_zip(genome_obj, project_file_cache)

    missing = len(
        [x for x in genome_status.values() if len(x.filename) == 0])
    logger.info(
        'Dataset has {} missing sets of antiSMASH data (from a total of {})'
        .format(missing, len(genome_records)))

    with open(genome_status_file, 'w', newline='\n') as f:
        for obj in genome_status.values():
            f.write(obj.to_csv() + '\n')

    if missing == len(genome_records):
        logger.warning('Failed to successfully retrieve ANY genome data!')
