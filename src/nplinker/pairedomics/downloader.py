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

import json
import os
import sys
import shutil
import zipfile
from deprecated import deprecated
import httpx
from progress.spinner import Spinner
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.gnps.gnps_downloader import GNPSDownloader
from nplinker.metabolomics.gnps.gnps_extractor import GNPSExtractor
from nplinker.strains import Strain
from nplinker.strain_collection import StrainCollection
from nplinker.genomics.mibig import download_and_extract_mibig_metadata
from nplinker.genomics.antismash import download_antismash_data


logger = LogConfig.getLogger(__name__)

from .runbigscape import run_bigscape


PAIREDOMICS_PROJECT_DATA_ENDPOINT = 'https://pairedomicsdata.bioinformatics.nl/api/projects'
PAIREDOMICS_PROJECT_URL = 'https://pairedomicsdata.bioinformatics.nl/api/projects/{}'
GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'

MIBIG_METADATA_URL = 'https://dl.secondarymetabolites.org/mibig/mibig_json_{}.tar.gz'
MIBIG_BGC_METADATA_URL = 'https://mibig.secondarymetabolites.org/repository/{}/annotations.json'
# MIBIG_BGC_GENBANK_URL = 'https://mibig.secondarymetabolites.org/repository/{}/{}.gbk'
# MIBIG_BGC_JSON_URL = 'https://mibig.secondarymetabolites.org/repository/{}/{}.json'


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

        self._init_folder_structure(local_cache)

        # init project json files
        self.all_project_json = None
        if not os.path.exists(self.project_json_file) or force_download:
            logger.info('Downloading new copy of platform project data...')
            self.all_project_json = self._download_and_load_json(
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
        self.project_json = self._download_and_load_json(
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
    
    def _init_folder_structure(self, local_cache):
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

        self._download_metabolomics_zipfile(self.gnps_task_id)

        download_antismash_data(self.project_json['genomes'])

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

    def _download_metabolomics_zipfile(self, gnps_task_id):
        archive = GNPSDownloader(gnps_task_id, self.project_download_cache).download().get_download_path()
        GNPSExtractor(archive, self.project_file_cache).extract()


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


    def _download_and_load_json(self, url, local_path):
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
