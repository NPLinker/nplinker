import sys
import os
import zipfile
import json

import httpx
from bs4 import BeautifulSoup
from xdg import XDG_CONFIG_HOME
from progress.bar import Bar
from progress.spinner import Spinner

from ..logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class Downloader(object):

    PAIREDOMICS_PROJECT_DATA_ENDPOINT = 'http://pairedomicsdata.bioinformatics.nl/api/projects'
    GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'
    ANTISMASH_DB_PAGE_URL = 'https://antismash-db.secondarymetabolites.org/output/{}'
    ANTISMASH_DB_DOWNLOAD_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/{}'

    NCBI_GENBANK_LOOKUP_URL = 'https://www.ncbi.nlm.nih.gov/nuccore/{}?report=docsum'
    NCBI_ASSEMBLY_LOOKUP_URL = 'https://www.ncbi.nlm.nih.gov/assembly?LinkName=nuccore_assembly&from_uid={}'

    def __init__(self, project_id, force_download=False):
        self.gnps_massive_id = project_id
        self.pairedomics_id = None
        self.gnps_task_id = None
        self.local_cache = os.path.join(os.getenv('HOME'), 'nplinker_data', 'pairedomics')
        self.local_download_cache = os.path.join(self.local_cache, 'downloads')
        self.local_file_cache = os.path.join(self.local_cache, 'extracted')
        self.local_project_json_file = os.path.join(self.local_cache, 'pairedomics_project_data.json')
        self.local_project_json_data = None
        os.makedirs(self.local_cache, exist_ok=True)

        self.json_data = None
        
        logger.info('Downloader for {}, caching to {}'.format(project_id, self.local_cache))

        if not os.path.exists(self.local_project_json_file) or force_download:
            logger.info('Downloading new copy of platform project data...')
            self._download_platform_project_data(self.local_project_json_file)
        else:
            logger.info('Using existing copy of platform project data')
            with open(self.local_project_json_file, 'r') as f:
                self.local_project_json_data = json.load(f)

        # find the specified project
        for project in self.local_project_json_data['data']:
            pairedomics_id = project['_id']
            gnps_massive_id = project['project']['metabolomics']['project']['GNPSMassIVE_ID']

            if gnps_massive_id == project_id:
                self.pairedomics_id = pairedomics_id
                if 'molecular_network' not in project['project']['metabolomics']['project']:
                    raise Exception('Dataset has no GNPS data URL!')

                self.gnps_task_id = project['project']['metabolomics']['project']['molecular_network']
                logger.info('Found project with internal ID {}'.format(self.pairedomics_id))
                self.project_json = project

        if self.pairedomics_id is None:
            raise Exception('Failed to find pairedomics project with ID {}'.format(self.gnps_massive_id))

        # create local cache folders for this dataset
        self.project_download_cache = os.path.join(self.local_download_cache, self.gnps_massive_id)
        os.makedirs(self.project_download_cache, exist_ok=True)

        self.project_file_cache = os.path.join(self.local_file_cache, self.gnps_massive_id)
        os.makedirs(self.project_file_cache, exist_ok=True)

    def get(self):
        logger.info('Going to download the metabolomics data file')

        self._download_metabolomics_zipfile(self.gnps_task_id)
        self._parse_genome_labels(self.project_json['project']['genome_metabolome_links'])
        self._download_genomics_data(self.project_json['project']['genomes'])

    def _is_new_gnps_format(self, directory):
        return os.path.exists(os.path.join(directory, 'qiime2_output'))

    def _download_genomics_data(self, genome_records):
        found = 0
        
        missing_cache_file = os.path.join(self.project_download_cache, 'missing_antismash.txt')
        missing_files = set()
        if os.path.exists(missing_cache_file):
            with open(missing_cache_file, 'r') as f:
                for l in f.readlines():
                    missing_files.add(l.strip())

        logger.debug('Dataset has {} missing sets of antiSMASH data'.format(len(missing_files)))

        for i, gr in enumerate(genome_records):
            label = gr['genome_label']
            if 'RefSeq_accession' in gr['genome_ID']:
                accession = gr['genome_ID']['RefSeq_accession']
                # TODO does original value need preserved
                accession = accession[:accession.rindex('.')]
                if accession in missing_files:
                    logger.warning('Not attempting to download data for accession={}'.format(accession))
                    continue
            
                logger.info('Checking for antismash data {}/{}, accession={}'.format(i+1, len(genome_records), accession))
                if self._download_antismash_zip(accession):
                    found += 1
                else:
                    missing_files.add(accession)

            elif 'GenBank_accession' in gr['genome_ID']:
                accession = gr['genome_ID']['GenBank_accession']
                try:
                    refseq_accession = self._get_refseq_from_genbank(accession)
                except:
                    logger.warning('Failed resolving GenBank accession {}'.format(accession))
                    continue

                if refseq_accession in missing_files:
                    logger.warning('Not attempting to download data for accession={}'.format(refseq_accession))
                    continue

                logger.info('Checking for antismash data {}/{}, accession={}'.format(i+1, len(genome_records), refseq_accession))
                if self._download_antismash_zip(refseq_accession):
                    found += 1
                else:
                    missing_files.add(refseq_accession)

            else:
                logger.warning('Missing RefSeq_accession label for {}'.format(label))
                continue

        with open(missing_cache_file, 'w') as f:
            for mf in missing_files:
                f.write('{}\n'.format(mf))

        print('Obtained {}/{} sets of genome data'.format(found, len(genome_records)))
        if found == 0:
            raise Exception('Failed to download ANY genome data!')

    def _get_refseq_from_genbank(self, genbank_id):
        """
        Super hacky way of trying to resolve the GenBank accession into RefSeq accession.
        """
        logger.info('Attempting to resolve RefSeq accession from Genbank accession {}'.format(genbank_id))
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

        # Look up genbank ID
        url = Downloader.NCBI_GENBANK_LOOKUP_URL.format(genbank_id)
        resp = httpx.get(url)
        soup = BeautifulSoup(resp.content, 'html.parser')
        ids = soup.find('dl', {'class': 'rprtid'})
        for field_idx, field in enumerate(ids.findChildren()):
            if field.getText().strip() == 'GI:':
                seq_id = ids.findChildren()[field_idx + 1].getText().strip()
                break

        # Look up assembly
        url = Downloader.NCBI_ASSEMBLY_LOOKUP_URL.format(seq_id)
        resp = httpx.get(url)
        soup = BeautifulSoup(resp.content, 'html.parser')
        title_href = soup.find('p', {'class': 'title'}).a['href']
        refseq_id = title_href.split('/')[-1].split('.')[0]

        return refseq_id


    def _download_antismash_zip(self, accession_id):
        # save files as <accession>.zip to avoid having to repeat above lookup every time
        local_path = os.path.join(self.project_download_cache, '{}.zip'.format(accession_id))
        logger.debug('Checking for existing antismash zip at {}'.format(local_path))
        cached = False
        if os.path.exists(local_path):
            logger.info('Found cached file at {}'.format(local_path))
            try:
                azip = zipfile.ZipFile(local_path)
                cached = True
            except zipfile.BadZipFile as bzf:
                logger.info('Invalid antismash zipfile found, will download again')
                os.unlink(local_path)

        if not cached:
            url = Downloader.ANTISMASH_DB_PAGE_URL.format(accession_id)
            logger.info('antismash DB lookup for {} => {}'.format(accession_id, url))
            resp = httpx.get(url)
            soup = BeautifulSoup(resp.content, 'html.parser')
            dl_list = soup.find('ul', {'id': 'downloadoptions'})
            if dl_list is None:
                logger.warning('Failed to download antismash-db results for {}'.format(accession_id))
                return False

            filename = dl_list.li.a['href']
            zipfile_url = Downloader.ANTISMASH_DB_DOWNLOAD_URL.format(accession_id, filename)
            with open(local_path, 'wb') as f:
                total_bytes, last_total = 0, 0
                with httpx.stream('GET', zipfile_url) as r:
                    logger.debug('zipfile URL is {}'.format(zipfile_url))
                    filesize = int(r.headers['content-length'])
                    bar = Bar(filename, max=filesize, suffix='%(percent)d%%')
                    for data in r.iter_bytes():
                        f.write(data)
                        total_bytes += len(data)
                        bar.next(len(data))
                    bar.finish()
        
        logger.debug('Extracting antismash data')

        antismash_zip = zipfile.ZipFile(local_path)
        kc_prefix = '{}/knownclusterblast'.format(accession_id)
        for zip_member in antismash_zip.namelist():
            # TODO other files here?
            if zip_member.endswith('.gbk') or zip_member.startswith(kc_prefix):
                antismash_zip.extract(zip_member, path=os.path.join(self.project_file_cache, 'antismash'))

        return True

    def _parse_genome_labels(self, genome_metabolome_links):
        # TODO have to update this after actually parsing metabolomics TSV
        for gml in genome_metabolome_links:
            label = gml['genome_label']
            filename = os.path.split(gml['metabolomics_file'])[1]

    def _download_metabolomics_zipfile(self, gnps_task_id):
        url = Downloader.GNPS_DATA_DOWNLOAD_URL.format(gnps_task_id)

        self.metabolomics_zip = os.path.join(self.project_download_cache, 'metabolomics_data.zip')

        cached = False
        if os.path.exists(self.metabolomics_zip):
            logger.info('Found existing metabolomics_zip at {}'.format(self.metabolomics_zip))
            try:
                mbzip = zipfile.ZipFile(self.metabolomics_zip)
                cached = True
            except zipfile.BadZipFile as bzf:
                logger.info('Invalid metabolomics zipfile found, will download again!')
                os.unlink(self.metabolomics_zip)
        
        if not cached:
            logger.info('Downloading metabolomics data from {}'.format(url))
            with open(self.metabolomics_zip, 'wb') as f:
                # note that this requires a POST, not a GET
                total_bytes, last_total = 0, 0
                spinner = Spinner('Downloading metabolomics data... ')
                with httpx.stream('POST', url) as r:
                    for data in r.iter_bytes():
                        f.write(data)
                        total_bytes += len(data)
                        spinner.next()
                spinner.finish()

        logger.info('Downloaded metabolomics data!')

        # this should throw an exception if zip is malformed etc
        mbzip = zipfile.ZipFile(self.metabolomics_zip)

        logger.info('Extracting files to {}'.format(self.project_file_cache))
        # extract the contents to the file cache folder. only want some of the files
        # so pick them out and only extract those:
        # - .mgf in root
        # - root/clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.tsv
        # - root/networkedges_selfloop/*.pairsinfo
        for member in mbzip.namelist():
            if member.endswith('.mgf')\
                or member.startswith('clusterinfosummarygroup_attributes_withIDs_withcomponentID')\
                or member.startswith('networkedges_selfloop'):
                    mbzip.extract(member, path=self.project_file_cache)
        if self._is_new_gnps_format(self.project_file_cache):
            logger.info('Found NEW GNPS structure')
        else:
            logger.info('Found OLD GNPS structure')

    def _download_platform_project_data(self, local_path):
        resp = httpx.get(Downloader.PAIREDOMICS_PROJECT_DATA_ENDPOINT)
        if not resp.status_code == 200:
            raise Exception('Failed to download pairedomics project JSON data! {}'.format(resp.status))

        logger.debug('Downloaded {} bytes of project data'.format(len(resp.content)))
        self.local_project_json_data = json.loads(resp.content)
        logger.info('Number of projects found: {}'.format(len(self.local_project_json_data['data'])))

        with open(local_path, 'w') as f:
           json.dump(self.local_project_json_data, f) 

        logger.info('Stored data in {}'.format(local_path))
    

if __name__ == "__main__":
    # salinispora dataset (only one that actually works)
    # d = Downloader('MSV000079284')
    d = Downloader('MSV000084771')
    d.get()

    # all_ids = ['MSV000078995',
    #         'MSV000082831',
    #         'MSV000078847',
    #         'MSV000078850',
    #         'MSV000084723',
    #         'MSV000084781',
    #         'MSV000078836',
    #         'MSV000078839',
    #         'MSV000084278',
    #         'MSV000079284',
    #         'MSV000084771']
    # failed = []
    # for id in all_ids:
    #     print('Attempting to get dataset {}'.format(id))
    #     try:
    #         d = Downloader(id)
    #         d.get()
    #     except Exception as e:
    #         print('Failed: {}'.format(e))
    #         failed.append(id)
    #         continue

    # print('Failed IDs')
    # print(failed)
