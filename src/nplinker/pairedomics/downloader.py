import json
import os
import shutil
import sys
import httpx
from nplinker.genomics.mibig import download_and_extract_mibig_metadata
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.gnps.gnps_downloader import GNPSDownloader
from nplinker.metabolomics.gnps.gnps_extractor import GNPSExtractor
from . import podp_download_and_extract_antismash_data
from .runbigscape import podp_run_bigscape


logger = LogConfig.getLogger(__name__)

PAIREDOMICS_PROJECT_DATA_ENDPOINT = 'https://pairedomicsdata.bioinformatics.nl/api/projects'
PAIREDOMICS_PROJECT_URL = 'https://pairedomicsdata.bioinformatics.nl/api/projects/{}'
GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'

MIBIG_METADATA_URL = 'https://dl.secondarymetabolites.org/mibig/mibig_json_{}.tar.gz'
MIBIG_BGC_METADATA_URL = 'https://mibig.secondarymetabolites.org/repository/{}/annotations.json'


class PODPDownloader():
    # TODO: move to independent config file  ---C.Geng
    PFAM_PATH = os.path.join(sys.prefix, 'nplinker_lib')

    def __init__(self, platform_id, force_download=False, local_cache=None):
        self.gnps_massive_id = platform_id
        self.pairedomics_id = None
        self.gnps_task_id = None
        self.json_data = None

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
            with open(self.all_project_json_file, encoding="utf-8") as f:
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
                logger.debug('platform_id %s matched to pairedomics_id %s',
                             self.gnps_massive_id, self.pairedomics_id)
                break

        if self.pairedomics_id is None:
            raise Exception(
                f'Failed to find a pairedomics project with ID {self.gnps_massive_id}'
            )

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

        with open(os.path.join(self.project_file_cache, 'platform_data.json'),
                  'w',
                  encoding='utf-8') as f:
            f.write(str(self.project_json))

    def _init_folder_structure(self, local_cache):
        """Create local cache folders and set up paths for various files"""

        # init local cache root
        self.local_cache = local_cache
        self.local_download_cache = os.path.join(self.local_cache, 'downloads')
        self.local_file_cache = os.path.join(self.local_cache, 'extracted')
        os.makedirs(self.local_cache, exist_ok=True)
        logger.info('PODPDownloader for %s, caching to %s',
                    self.gnps_massive_id, self.local_cache)

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

        # init project paths
        self.all_project_json_file = os.path.join(self.local_cache,
                                                  'all_projects.json')
        self.project_json_file = os.path.join(self.local_cache,
                                              f'{self.gnps_massive_id}.json')

    # download function
    def get(self, do_bigscape, extra_bigscape_parameters, use_mibig,
            mibig_version):
        logger.info('Going to download the metabolomics data file')

        self._download_metabolomics_zipfile(self.gnps_task_id)

        # TODO CG: this function will modify the project_json['genomes'],
        # this should be done in a better way
        podp_download_and_extract_antismash_data(self.project_json['genomes'],
                                                 self.project_download_cache,
                                                 self.project_file_cache)

        if use_mibig:
            self._download_mibig_json(mibig_version)
        podp_run_bigscape(self.project_file_cache, self.PFAM_PATH, do_bigscape,
                          extra_bigscape_parameters)

    def _download_mibig_json(self, version):
        output_path = os.path.join(self.project_file_cache, 'mibig_json')

        # Override existing mibig json files
        if os.path.exists(output_path):
            shutil.rmtree(output_path)

        os.makedirs(output_path)

        download_and_extract_mibig_metadata(self.project_download_cache,
                                            output_path, version)

        self._create_completed_file(output_path)

        return True

    @staticmethod
    def _create_completed_file(output_path):
        with open(os.path.join(output_path, 'completed'),
                  'w',
                  encoding='utf-8'):
            pass


    def _download_metabolomics_zipfile(self, gnps_task_id):
        archive = GNPSDownloader(
            gnps_task_id,
            self.project_download_cache).download().get_download_path()
        GNPSExtractor(archive, self.project_file_cache).extract()

    def _download_and_load_json(self, url, local_path):
        resp = httpx.get(url, follow_redirects=True)
        if not resp.status_code == 200:
            raise Exception(
                f'Failed to download {url} (status code {resp.status_code})')

        content = json.loads(resp.content)
        with open(local_path, 'w', encoding='utf-8') as f:
            json.dump(content, f)

        logger.debug('Downloaded %s to %s', url, local_path)

        return content
