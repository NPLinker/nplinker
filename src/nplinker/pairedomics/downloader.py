import json
import os
from os import PathLike
from pathlib import Path
import shutil
from nplinker.genomics.mibig import download_and_extract_mibig_metadata
from nplinker.globals import PFAM_PATH
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.gnps.gnps_downloader import GNPSDownloader
from nplinker.metabolomics.gnps.gnps_extractor import GNPSExtractor
from nplinker.utils import download_url
from . import podp_download_and_extract_antismash_data
from .runbigscape import podp_run_bigscape


logger = LogConfig.getLogger(__name__)

PAIREDOMICS_PROJECT_DATA_ENDPOINT = 'https://pairedomicsdata.bioinformatics.nl/api/projects'
PAIREDOMICS_PROJECT_URL = 'https://pairedomicsdata.bioinformatics.nl/api/projects/{}'
GNPS_DATA_DOWNLOAD_URL = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task={}&view=download_clustered_spectra'

MIBIG_METADATA_URL = 'https://dl.secondarymetabolites.org/mibig/mibig_json_{}.tar.gz'
MIBIG_BGC_METADATA_URL = 'https://mibig.secondarymetabolites.org/repository/{}/annotations.json'


class PODPDownloader():

    def __init__(self,
                 podp_platform_id: str,
                 force_download: bool = False,
                 root_dir: str | PathLike | None = None):
        """Downloader for PODP pipeline.

        The downloader will download the following data:
            - GNPS Molecular Network task results
            - AntiSMASH results
            - MIBiG metadata

        Args:
            podp_platform_id(str): The metabolomics project ID of PODP platform,
                e.g. GNPS MassIVE ID.
            force_download (bool): Re-download data even if it already exists
                locally. Defaults to False.
            working_dir (str | PathLike | None): The root directory to use for
                the project. Defaults to None, in which case the default location
                is used.

        Raises:
            ValueError: If the given ID does not have a corresponding PODP ID,
                or if the GNPS Molecular Network task URL does not exist for
                the given ID.
        """
        self.gnps_massive_id = podp_platform_id

        if root_dir is None:
            root_dir = os.path.join(os.getenv('HOME'), 'nplinker_data',
                                    'pairedomics')

        # TODO CG: init folder structure should be moved out of PODPDownloader
        self._init_folder_structure(root_dir)

        # init project json files
        if not os.path.exists(self.project_json_file) or force_download:
            logger.info('Downloading new copy of platform project data...')
            self.all_projects_json_data = self._download_and_load_json(
                PAIREDOMICS_PROJECT_DATA_ENDPOINT, self.all_projects_json_file)
        else:
            logger.info('Using existing copy of platform project data')
            with open(self.all_projects_json_file, encoding="utf-8") as f:
                self.all_projects_json_data = json.load(f)

        # Verify that the given ID has a corresponding PODP ID
        self.podp_id = None
        for project in self.all_projects_json_data['data']:
            if self.gnps_massive_id == project['metabolite_id']:
                self.podp_id = project['_id']
                logger.debug('Given ID %s matched to PODP ID %s',
                             self.gnps_massive_id, self.podp_id)
                break
        if self.podp_id is None:
            raise ValueError(
                f'Failed to find PODP ID for given ID {self.gnps_massive_id}')

        # now get the project JSON data
        logger.info('Found project, retrieving JSON data...')
        self.project_json_data = self._download_and_load_json(
            PAIREDOMICS_PROJECT_URL.format(self.podp_id),
            self.project_json_file)

        self.gnps_task_id = self.project_json_data['metabolomics'][
            'project'].get('molecular_network')
        if self.gnps_task_id is None:
            raise ValueError(
                f'GNPS Molecular Network task URL not exist for '
                f'given ID {self.gnps_massive_id}. Please check and'
                f'run GNPS Molecular Network task first.')

    def _init_folder_structure(self, working_dir):
        """Create local cache folders and set up paths for various files"""

        # init local cache root
        self.working_dir = working_dir
        self.downloads_dir = os.path.join(self.working_dir, 'downloads')
        self.results_dir = os.path.join(self.working_dir, 'extracted')
        os.makedirs(self.working_dir, exist_ok=True)
        logger.info('PODPDownloader for %s, caching to %s',
                    self.gnps_massive_id, self.working_dir)

        # create local cache folders for this dataset
        self.project_downloads_dir = os.path.join(self.downloads_dir,
                                                  self.gnps_massive_id)
        os.makedirs(self.project_downloads_dir, exist_ok=True)

        self.project_results_dir = os.path.join(self.results_dir,
                                                self.gnps_massive_id)
        os.makedirs(self.project_results_dir, exist_ok=True)

        # placeholder directories
        for d in ['antismash', 'bigscape']:
            os.makedirs(os.path.join(self.project_results_dir, d),
                        exist_ok=True)

        # init project paths
        self.all_projects_json_file = os.path.join(self.working_dir,
                                                   'all_projects.json')
        self.project_json_file = os.path.join(self.working_dir,
                                              f'{self.gnps_massive_id}.json')

    # download function
    def get(self, do_bigscape, extra_bigscape_parameters, use_mibig,
            mibig_version):
        logger.info('Going to download the metabolomics data file')

        self._download_metabolomics_zipfile(self.gnps_task_id)

        # TODO CG: this function will modify the project_json['genomes'],
        # this should be done in a better way
        podp_download_and_extract_antismash_data(
            self.project_json_data['genomes'], self.project_downloads_dir,
            self.project_results_dir)

        if use_mibig:
            self._download_mibig_json(mibig_version)
        podp_run_bigscape(self.project_results_dir, PFAM_PATH, do_bigscape,
                          extra_bigscape_parameters)

    def _download_mibig_json(self, version):
        output_path = os.path.join(self.project_results_dir, 'mibig_json')

        # Override existing mibig json files
        if os.path.exists(output_path):
            shutil.rmtree(output_path)

        os.makedirs(output_path)

        download_and_extract_mibig_metadata(self.project_downloads_dir,
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
            self.project_downloads_dir).download().get_download_file()
        GNPSExtractor(archive, self.project_results_dir).extract()

    def _download_and_load_json(self, url: str,
                                output_file: str | PathLike) -> dict:
        """Download a JSON file from a URL and return the parsed JSON data"""
        fpath = Path(output_file)
        download_url(url, fpath.parent, fpath.name)
        logger.debug('Downloaded %s to %s', url, output_file)

        with open(output_file, 'r') as f:
            data = json.load(f)

        return data
