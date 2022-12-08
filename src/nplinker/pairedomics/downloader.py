# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import csv
import io
import json
import os
from pathlib import Path
import re
import sys
import shutil
import tarfile
import time
import zipfile
from deprecated import deprecated
import httpx
from bs4 import BeautifulSoup
from progress.bar import Bar
from progress.spinner import Spinner
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.gnps.gnps_downloader import GNPSDownloader
from nplinker.metabolomics.gnps.gnps_extractor import GNPSExtractor
from nplinker.strains import Strain
from nplinker.strain_collection import StrainCollection
from nplinker.genomics.mibig import download_and_extract_mibig_metadata


logger = LogConfig.getLogger(__name__)

from .runbigscape import run_bigscape


PAIREDOMICS_PROJECT_DATA_ENDPOINT = 'https://pairedomicsdata.bioinformatics.nl/api/projects'
PAIREDOMICS_PROJECT_URL = 'https://pairedomicsdata.bioinformatics.nl/api/projects/{}'
GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'

ANTISMASH_DB_PAGE_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/'
ANTISMASH_DB_DOWNLOAD_URL = 'https://antismash-db.secondarymetabolites.org/output/{}/{}'

ANTISMASH_DBV2_PAGE_URL = 'https://antismash-dbv2.secondarymetabolites.org/output/{}/'
ANTISMASH_DBV2_DOWNLOAD_URL = 'https://antismash-dbv2.secondarymetabolites.org/output/{}/{}'

NCBI_LOOKUP_URL_NEW = 'https://www.ncbi.nlm.nih.gov/assembly/?term={}'

JGI_GENOME_LOOKUP_URL = 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid={}'

MIBIG_METADATA_URL = 'https://dl.secondarymetabolites.org/mibig/mibig_json_{}.tar.gz'
MIBIG_BGC_METADATA_URL = 'https://mibig.secondarymetabolites.org/repository/{}/annotations.json'
# MIBIG_BGC_GENBANK_URL = 'https://mibig.secondarymetabolites.org/repository/{}/{}.gbk'
# MIBIG_BGC_JSON_URL = 'https://mibig.secondarymetabolites.org/repository/{}/{}.json'

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


class Downloader():
    # TODO: move to independent config file  ---C.Geng
    PFAM_PATH = os.path.join(sys.prefix, 'nplinker_lib')

    def __init__(self, platform_id, force_download=False, local_cache = None):
        self.gnps_massive_id = platform_id
        self.pairedomics_id = None
        self.gnps_task_id = None
        self.json_data = None        
        self.strains = StrainCollection()
        self.growth_media = {}

        if local_cache is None:
            local_cache = os.path.join(os.getenv('HOME'), 'nplinker_data',
                                        'pairedomics')

        self.init_folder_structure(local_cache)

        # init project json files
        self.all_project_json = None
        if not os.path.exists(self.project_json_file) or force_download:
            logger.info('Downloading new copy of platform project data...')
            self.all_project_json = self._download_platform_json_to_file(
                PAIREDOMICS_PROJECT_DATA_ENDPOINT, self.all_project_json_file)
        else:
            logger.info('Using existing copy of platform project data')
            with open(self.all_project_json_file) as f:
                self.all_project_json = json.load(f)

        # query the pairedomics webservice with the project ID to retrieve the data. unfortunately
        # this is not the MSV... ID, but an internal GUID string. To get that, first need to get the
        # list of all projects, find the one with a 'metabolite_id' value matching the MSV... ID, and
        # then extract its '_id' value to get the GUID

        # find the specified project and store its ID
        for project in self.all_project_json['data']:
            pairedomics_id = project['_id']
            gnps_massive_id = project['metabolite_id']

            if gnps_massive_id == platform_id:
                self.pairedomics_id = pairedomics_id
                logger.debug(
                    'platform_id {} matched to pairedomics_id {}'.format(
                        self.gnps_massive_id, self.pairedomics_id))
                break

        if self.pairedomics_id is None:
            raise Exception(
                'Failed to find a pairedomics project with ID {}'.format(
                    self.gnps_massive_id))

        # now get the project JSON data
        self.project_json = None
        logger.info('Found project, retrieving JSON data...')
        self.project_json = self._download_platform_json_to_file(
            PAIREDOMICS_PROJECT_URL.format(self.pairedomics_id),
            self.project_json_file)

        if 'molecular_network' not in self.project_json['metabolomics'][
                'project']:
            raise Exception('Dataset has no GNPS data URL!')

        self.gnps_task_id = self.project_json['metabolomics']['project'][
            'molecular_network']

        with open(os.path.join(self.project_file_cache,
                                  'platform_data.json'),
                     'w',
                     encoding='utf-8') as f:
            f.write(str(self.project_json))
    
    def init_folder_structure(self, local_cache):
        # init local cache root   
        self.local_cache = local_cache
        self.local_download_cache = os.path.join(self.local_cache, 'downloads')
        self.local_file_cache = os.path.join(self.local_cache, 'extracted')    
        os.makedirs(self.local_cache, exist_ok=True)
        logger.info('Downloader for {}, caching to {}'.format(
            self.gnps_massive_id, self.local_cache))

        # create local cache folders for this dataset
        self.project_download_cache = os.path.join(self.local_download_cache,
                                                   self.gnps_massive_id)
        os.makedirs(self.project_download_cache, exist_ok=True)

        self.project_file_cache = os.path.join(self.local_file_cache,
                                               self.gnps_massive_id)
        os.makedirs(self.project_file_cache, exist_ok=True)

        # placeholder directories
        for d in ['antismash', 'bigscape']:
            os.makedirs(os.path.join(self.project_file_cache, d),
                        exist_ok=True)
        
        # init strain mapping filepath
        self.strain_mappings_file = os.path.join(self.project_file_cache,
                                                 'strain_mappings.csv') 
        
        # init project paths
        self.all_project_json_file = os.path.join(self.local_cache,
                                                  'all_projects.json')
        self.project_json_file = os.path.join(
            self.local_cache, f'{self.gnps_massive_id}.json')


	# CG: download function
    def get(self, do_bigscape, extra_bigscape_parameters, use_mibig,
            mibig_version):
        logger.info('Going to download the metabolomics data file')

        self._download_metabolomics_zipfile_v2(self.gnps_task_id)

        self._download_genomics_data(self.project_json['genomes'])

        # CG: it extracts strain names and later will be used for strains
        self._parse_genome_labels(self.project_json['genome_metabolome_links'],
                                  self.project_json['genomes'])

        # CG: it generates the strain_mapping.csv file
        self.strains.generate_strain_mappings(self.strain_mappings_file,
            os.path.join(self.project_file_cache, 'antismash'))

        if use_mibig:
            self._download_mibig_json(mibig_version)
        self._run_bigscape(do_bigscape, extra_bigscape_parameters)

    def _is_new_gnps_format(self, directory):
        # TODO this should test for existence of quantification table instead
        return os.path.exists(os.path.join(directory, 'qiime2_output'))

    def _run_bigscape(self, do_bigscape, extra_bigscape_parameters):
        # TODO this currently assumes docker environment, allow customisation?
        # can check if in container with: https://stackoverflow.com/questions/20010199/how-to-determine-if-a-process-runs-inside-lxc-docker
        if not do_bigscape:
            logger.info('BiG-SCAPE disabled by configuration, not running it')
            return

        logger.info('Running BiG-SCAPE! extra_bigscape_parameters="{}"'.format(
            extra_bigscape_parameters))
        try:
            run_bigscape('bigscape.py',
                         os.path.join(self.project_file_cache, 'antismash'),
                         os.path.join(self.project_file_cache, 'bigscape'),
                         self.PFAM_PATH, extra_bigscape_parameters)
        except Exception as e:
            logger.warning(
                'Failed to run BiG-SCAPE on antismash data, error was "{}"'.
                format(e))

    def _ncbi_genbank_search(self, genbank_id, retry_time=5.0):
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

    def _resolve_genbank_accession(self, genbank_id):
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
            dl_element = self._ncbi_genbank_search(genbank_id)
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

    def _resolve_jgi_accession(self, jgi_id):
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

        return self._resolve_genbank_accession(link.text)

    def _get_best_available_genome_id(self, genome_id_data):
        if 'RefSeq_accession' in genome_id_data:
            return genome_id_data['RefSeq_accession']
        elif 'GenBank_accession' in genome_id_data:
            return genome_id_data['GenBank_accession']
        elif 'JGI_Genome_ID' in genome_id_data:
            return genome_id_data['JGI_Genome_ID']

        logger.warning('No known genome ID field in genome data: {}'.format(
            genome_id_data))
        return None

    def _resolve_genome_id_data(self, genome_id_data):
        if 'RefSeq_accession' in genome_id_data:
            # best case, can use this directly
            return genome_id_data['RefSeq_accession']
        elif 'GenBank_accession' in genome_id_data:
            # resolve via NCBI
            return self._resolve_genbank_accession(
                genome_id_data['GenBank_accession'])
        elif 'JGI_Genome_ID' in genome_id_data:
            # resolve via JGI => NCBI
            return self._resolve_jgi_accession(genome_id_data['JGI_Genome_ID'])

        logger.warning(
            f'Unable to resolve genome_ID: {genome_id_data}')
        return None

    def _download_genomics_data(self, genome_records):
        genome_status = {}

        # this file records genome IDs and local filenames to avoid having to repeat HTTP requests
        # each time the app is loaded (this can take a lot of time if there are dozens of genomes)
        genome_status_file = os.path.join(self.project_download_cache,
                                          'genome_status.txt')

        # genome lookup status info
        if os.path.exists(genome_status_file):
            with open(genome_status_file) as f:
                for line in csv.reader(f):
                    asobj = GenomeStatus.from_csv(*line)
                    genome_status[asobj.original_id] = asobj

        for i, genome_record in enumerate(genome_records):
            label = genome_record['genome_label']

            # get the best available ID from the dict
            best_id = self._get_best_available_genome_id(
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

                genome_obj.resolved_id = self._resolve_genome_id_data(
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
                if self._download_antismash_zip(genome_obj):
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

            self._extract_antismash_zip(genome_obj)

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

    def _download_mibig_json(self, version):
        output_path = os.path.join(self.project_file_cache, 'mibig_json')

        # Override existing mibig json files
        if os.path.exists(output_path):
            shutil.rmtree(output_path)

        os.makedirs(output_path)

        download_and_extract_mibig_metadata(self.project_download_cache,
                                        output_path, version)

        open(os.path.join(output_path, 'completed'), 'w').close()

        return True

    def _get_antismash_db_page(self, genome_obj):
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

    def _get_antismash_zip_data(self, accession_id, filename, local_path):
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

    def _download_antismash_zip(self, antismash_obj):
        # save zip files to avoid having to repeat above lookup every time
        local_path = os.path.join(self.project_download_cache,
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
            filename = self._get_antismash_db_page(antismash_obj)
            if filename is None:
                return False

            self._get_antismash_zip_data(antismash_obj.resolved_id, filename,
                                         local_path)
            antismash_obj.filename = local_path

        return True

    def _extract_antismash_zip(self, antismash_obj):
        if antismash_obj.filename is None or len(antismash_obj.filename) == 0:
            return False

        output_path = os.path.join(self.project_file_cache, 'antismash',
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

    def _parse_genome_labels(self, met_records, gen_records):
        temp = {}
        mc, gc = 0, 0

        # this method is supposed to extract the fields from the JSON data
        # which map strain names to mzXML files on the metabolomics side,
        # and to BGCs on the genomics side, and use that data to build a set
        # of NPLinker Strain objects for the current dataset

        # metabolomics: each of the JSON records should contain a field named
        # "genome_label", which is the one that should be used as the canonical
        # name for this strain by nplinker. Another field is called "metabolomics_file",
        # and this contains a URL to the corresponding mzXML file. so we want to
        # create a set of mappings from one to the other, with the complication that
        # there might be mappings from 2 or more mzXMLs to a single strain.
        # also should record the growth medium using the "sample_preparation_label" field.
        for rec in met_records:
            # this is the global strain identifier we should use
            label = rec['genome_label']
            # only want to record the actual filename of the mzXML URL
            filename = os.path.split(rec['metabolomics_file'])[1]

            # add the mzXML mapping for this strain
            if label in temp:
                temp[label].append(filename)
            else:
                temp[label] = [filename]
            mc += 1

            if label in self.growth_media:
                self.growth_media[label].add(rec['sample_preparation_label'])
            else:
                self.growth_media[label] = {
                    rec['sample_preparation_label']}

        for rec in gen_records:
            label = rec['genome_label']
            accession = rec.get('resolved_id', None)
            if accession is None:
                # this will happen for genomes where we couldn't retrieve data or resolve the ID
                logger.warning(
                    'Failed to extract accession from genome with label {}'.
                    format(label))
                continue

            if label in temp:
                temp[label].append(accession)
            else:
                temp[label] = [accession]
                gc += 1

        logger.info('Extracted {} strains from JSON (met={}, gen={})'.format(
            len(temp), mc, gc))
        for strain_label, strain_aliases in temp.items():
            strain = Strain(strain_label)
            for alias in strain_aliases:
                strain.add_alias(alias)
            self.strains.add(strain)

    def _download_metabolomics_zipfile_v2(self, gnps_task_id):
        archive = GNPSDownloader(gnps_task_id, self.project_download_cache).download().get_download_path()
        GNPSExtractor(archive, self.project_file_cache).extract()


    @deprecated
    def _download_metabolomics_zipfile(self, gnps_task_id):
        mbzip = self._load_gnps_data(gnps_task_id)
        self._extract_metabolomics_data(mbzip)
        self._log_gnps_format()


    @deprecated
    def _extract_metabolomics_data(self, mbzip):
        logger.info(f'Extracting files to {self.project_file_cache}')
        # extract the contents to the file cache folder. only want some of the files
        # so pick them out and only extract those:
        # - root/spectra/*.mgf
        # - root/clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.tsv
        # - root/networkedges_selfloop/*.pairsinfo
        # - root/quantification_table*
        # - root/metadata_table*
        # - root/DB_result*
        for member in mbzip.namelist():
            if member.startswith('clusterinfosummarygroup_attributes_withIDs_withcomponentID')\
                or member.startswith('networkedges_selfloop')\
                or member.startswith('quantification_table')\
                or member.startswith('metadata_table')\
                or member.startswith('DB_result')\
                or member.startswith('result_specnets_DB'):
                mbzip.extract(member, path=self.project_file_cache)
            # move the MGF file to a /spectra subdirectory to better fit expected structure
            elif member.endswith('.mgf'):
                os.makedirs(os.path.join(self.project_file_cache, 'spectra'),
                            exist_ok=True)
                mbzip.extract(member,
                              path=os.path.join(self.project_file_cache,
                                                'spectra'))

    @deprecated
    def _log_gnps_format(self):
        if self._is_new_gnps_format(self.project_file_cache):
            logger.info('Found NEW GNPS structure')
        else:
            logger.info('Found OLD GNPS structure')

    @deprecated
    def _load_gnps_data(self, gnps_task_id) -> zipfile.ZipFile:

        self.metabolomics_zip = os.path.join(self.project_download_cache,
                                             'metabolomics_data.zip')

        # Try read from cache
        if os.path.exists(self.metabolomics_zip):
            logger.info('Found existing metabolomics_zip at {}'.format(
                self.metabolomics_zip))
            try:
                mbzip = zipfile.ZipFile(self.metabolomics_zip)
                return mbzip
            except zipfile.BadZipFile as bzf:
                logger.info(
                    'Invalid metabolomics zipfile found, will download again!')
                os.unlink(self.metabolomics_zip)
        url = _generate_gnps_download_url(gnps_task_id)
        _execute_download(url, self.metabolomics_zip)

        # this should throw an exception if zip is malformed etc
        mbzip = zipfile.ZipFile(self.metabolomics_zip)
        return mbzip


    def _download_platform_json_to_file(self, url, local_path):
        resp = httpx.get(url, follow_redirects=True)
        if not resp.status_code == 200:
            raise Exception('Failed to download {} (status code {})'.format(
                url, resp.status_code))

        content = json.loads(resp.content)
        with open(local_path, 'w') as f:
            json.dump(content, f)

        logger.debug(f'Downloaded {url} to {local_path}')

        return content


@deprecated
def _generate_gnps_download_url(gnps_task_id):
    url = GNPS_DATA_DOWNLOAD_URL.format(gnps_task_id)
    return url

@deprecated
def _execute_download(url, metabolomics_zip):
    logger.info(f'Downloading metabolomics data from {url}')
    with open(metabolomics_zip, 'wb') as f:
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
