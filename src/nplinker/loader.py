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

import glob
import os
import sys
import time
import xml.etree.ElementTree as ET
from nplinker.annotations import load_annotations
from nplinker.class_info.chem_classes import ChemClassPredictions
from nplinker.class_info.class_matches import ClassMatches
from nplinker.class_info.runcanopus import run_canopus
from nplinker.genomics import loadBGC_from_cluster_files
from nplinker.genomics.mibig import MibigBGCLoader
from nplinker.logconfig import LogConfig
from nplinker.metabolomics import load_dataset
from nplinker.pairedomics.downloader import Downloader
from nplinker.genomics.mibig import download_and_extract_mibig_metadata
from nplinker.pairedomics.runbigscape import run_bigscape
from nplinker.strain_collection import StrainCollection


try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

logger = LogConfig.getLogger(__file__)


def find_via_glob(path, file_type, optional=False):
    try:
        filename = glob.glob(path)[0]
        return filename
    except (OSError, IndexError) as e:
        if not optional:
            # "from None" suppresses the traceback for the original exception, which isn't really needed
            raise Exception('ERROR: unable to find {} in path "{}"'.format(
                file_type, path)) from None

        logger.warn('WARNING: unable to find {} in path "{}"'.format(
            file_type, path))
        return None


def find_via_glob_alts_dir(paths, file_type, optional=False):
    path = None
    for p in paths:
        if os.path.exists(p):
            path = p
            break

    if path is None and not optional:
        raise Exception('ERROR: unable to find {} in {} paths: ({})'.format(
            file_type, len(paths), paths))
    elif path is None:
        logger.warning('WARNING: unable to find {} in {} paths: ({})'.format(
            file_type, len(paths), paths))

    return path


def find_via_glob_alts(paths, file_type, optional=False):
    filename = None
    for path in paths:
        try:
            filename = glob.glob(path)[0]
            break
        except (OSError, IndexError) as e:
            continue

    if filename is None and not optional:
        raise Exception('ERROR: unable to find {} in {} paths: ({})'.format(
            file_type, len(paths), paths))
    elif filename is None:
        logger.warning('WARNING: unable to find {} in {} paths: ({})'.format(
            file_type, len(paths), paths))

    return filename


def find_bigscape_dir(broot):
    logger.info(
        f'Trying to discover correct bigscape directory under {broot}')
    for root, dirs, files in os.walk(broot):
        if 'Network_Annotations_Full.tsv' in files:
            logger.info(f'Found network files directory: {root}')
            return root

    return None


class DatasetLoader():

    ANTISMASH_FMT_DEFAULT = 'default'
    ANTISMASH_FMT_FLAT = 'flat'
    ANTISMASH_FMTS = [ANTISMASH_FMT_DEFAULT, ANTISMASH_FMT_FLAT]
    ANTISMASH_DELIMITERS_DEFAULT = ['.', '_', '-']
    ANTISMASH_IGNORE_SPACES_DEFAULT = False

    TABLES_CUTOFF_DEFAULT = 2.0

    BIGSCAPE_CUTOFF_DEFAULT = 30
    EXTENDED_METADATA_TABLE_PARSING_DEFAULT = False
    USE_MIBIG_DEFAULT = True
    MIBIG_VERSION_DEFAULT = "1.4"

    RUN_BIGSCAPE_DEFAULT = True
    EXTRA_BIGSCAPE_PARAMS_DEFAULT = '--mibig --clans-off'

    RUN_CANOPUS_DEFAULT = False
    EXTRA_CANOPUS_PARAMS_DEFAULT = '--maxmz 600 formula zodiac structure canopus'

    # keys for overriding metabolomics data elements
    OR_NODES = 'nodes_file'
    OR_EDGES = 'edges_file'
    OR_EXTRA_NODES = 'extra_nodes_file'
    OR_MGF = 'mgf_file'
    OR_METADATA = 'metadata_table_file'
    OR_QUANT = 'quantification_table_file'
    OR_ANNO = 'annotations_dir'
    OR_ANNO_CONFIG = 'annotations_config_file'
    # and the same for genomics data
    OR_ANTISMASH = 'antismash_dir'
    OR_BIGSCAPE = 'bigscape_dir'
    OR_MIBIG_JSON = 'mibig_json_dir'
    OR_STRAINS = 'strain_mappings_file'
    # misc files
    OR_PARAMS = 'gnps_params_file'
    OR_DESCRIPTION = 'description_file'
    OR_INCLUDE_STRAINS = 'include_strains_file'
    # class predictions
    OR_CANOPUS = 'canopus_dir'
    OR_MOLNETENHANCER = 'molnetenhancer_dir'

    BIGSCAPE_PRODUCT_TYPES = [
        'NRPS', 'Others', 'PKSI', 'PKS-NRP_Hybrids', 'PKSother', 'RiPPs',
        'Saccharides', 'Terpene'
    ]

    # TODO: move to a config file, used by multiple modules
    PFAM_PATH = os.path.join(sys.prefix, 'nplinker_lib')

    def __init__(self, config_data):
        self._config = config_data
        self._dataset = config_data['dataset']
        self._docker = config_data.get('docker', {})
        self._webapp = config_data.get('webapp', {})
        self._antismash = config_data.get('antismash', {})
        self._overrides = self._dataset.get('overrides', {})
        self._antismash_delimiters = self._antismash.get(
            'antismash_delimiters', self.ANTISMASH_DELIMITERS_DEFAULT)
        self._antismash_format = self._antismash.get(
            'antismash_format', self.ANTISMASH_FMT_DEFAULT)
        self._antismash_ignore_spaces = self._antismash.get(
            'ignore_spaces', self.ANTISMASH_IGNORE_SPACES_DEFAULT)
        self._bigscape_cutoff = self._dataset.get('bigscape_cutoff',
                                                  self.BIGSCAPE_CUTOFF_DEFAULT)
        self._extended_metadata_table_parsing = self._dataset.get(
            'extended_metadata_table_parsing',
            self.EXTENDED_METADATA_TABLE_PARSING_DEFAULT)
        self._use_mibig = self._dataset.get('use_mibig',
                                            self.USE_MIBIG_DEFAULT)
        self._mibig_version = self._dataset.get('mibig_version',
                                                self.MIBIG_VERSION_DEFAULT)
        self._root = self._config['dataset']['root']
        self._platform_id = self._config['dataset']['platform_id']
        self._remote_loading = len(self._platform_id) > 0
        self.datadir = files('nplinker').joinpath('data')
        self.dataset_id = os.path.split(
            self._root)[-1] if not self._remote_loading else self._platform_id
        if self._remote_loading:
            self._downloader = Downloader(self._platform_id)
        else:
            self._downloader = None
        self.bgcs, self.gcfs, self.spectra, self.molfams = [], [], [], []
        self.mibig_bgc_dict = {}
        self.product_types = []
        logger.debug('DatasetLoader({}, {}, {})'.format(
            self._root, self._platform_id, self._remote_loading))

    def validate(self):
        # check antismash format is recognised
        if self._antismash_format not in self.ANTISMASH_FMTS:
            raise Exception('Unknown antismash format: {}'.format(
                self._antismash_format))

        # if remote loading mode, need to download the data here
        if self._remote_loading:
            self._root = self._downloader.project_file_cache
            logger.debug('remote loading mode, configuring root={}'.format(
                self._root))
            # CG: to download both MET and GEN data
            self._downloader.get(
                self._docker.get('run_bigscape', self.RUN_BIGSCAPE_DEFAULT),
                self._docker.get('extra_bigscape_parameters',
                                 self.EXTRA_BIGSCAPE_PARAMS_DEFAULT),
                self._use_mibig, self._mibig_version)

        # construct the paths and filenames required to load everything else and check
        # they all seem to exist (but don't parse anything yet)

        # 1. strain mapping are used for everything else so
        self.strain_mappings_file = self._overrides.get(
            self.OR_STRAINS) or os.path.join(self._root, 'strain_mappings.csv')

        # 2. MET: <root>/clusterinfo_summary/<some UID>.tsv (or .clustersummary apparently...) / nodes_file=<override>
        self.nodes_file = self._overrides.get(
            self.OR_NODES) or find_via_glob_alts([
                os.path.join(self._root, 'clusterinfo*', '*.tsv'),
                os.path.join(self._root, 'clusterinfo*', '*.clustersummary')
            ], self.OR_NODES)

        # 3. MET: <root>/networkedges_selfloop/<some UID>.selfloop (new) or .pairsinfo (old) / edges_file=<override>
        self.edges_file = self._overrides.get(
            self.OR_EDGES) or find_via_glob_alts([
                os.path.join(self._root, 'networkedges_selfloop',
                             '*.pairsinfo'),
                os.path.join(self._root, 'networkedges_selfloop', '*.selfloop')
            ], self.OR_EDGES)

        # 4. MET: <root>/*.csv / extra_nodes_file=<override>
        # TODO is the glob input OK?
        # => wait for updated dataset with latest output format
        # NOTE: only optional for Crusemann or Crusemann-like dataset format!
        self.extra_nodes_file = self._overrides.get(
            self.OR_EXTRA_NODES) or find_via_glob(os.path.join(
                self._root, 'quantification_table_reformatted', '*.csv'),
                                                  self.OR_EXTRA_NODES,
                                                  optional=True)

        # 5. MET: <root>/spectra/*.mgf (or <root>/*.mgf)/ mgf_file=<override>
        self.mgf_file = self._overrides.get(self.OR_MGF) or find_via_glob_alts(
            [
                os.path.join(self._root, 'spectra', '*.mgf'),
                os.path.join(self._root, '*.mgf')
            ], self.OR_MGF)

        # 6. MET: <root>/metadata_table/metadata_table-<number>.txt / metadata_table_file=<override>
        self.metadata_table_file = self._overrides.get(
            self.OR_METADATA) or find_via_glob(os.path.join(
                self._root, 'metadata_table', 'metadata_table*.txt'),
                                               self.OR_METADATA,
                                               optional=True)

        # 7. MET: <root>/quantification_table/quantification_table-<number>.csv / quantification_table_file=<override>
        self.quantification_table_file = self._overrides.get(
            self.OR_QUANT) or find_via_glob(os.path.join(
                self._root, 'quantification_table',
                'quantification_table*.csv'),
                                            self.OR_QUANT,
                                            optional=True)

        # 8. MET: <root>/DB_result/*.tsv (new) or <root>/result_specnets_DB/*.tsv (old) / annotations_dir=<override>
        self.annotations_dir = self._overrides.get(
            self.OR_ANNO) or find_via_glob_alts_dir([
                os.path.join(self._root, 'DB_result'),
                os.path.join(self._root, 'result_specnets_DB')
            ], self.OR_ANNO)
        self.annotations_config_file = self._overrides.get(
            self.OR_ANNO_CONFIG) or os.path.join(self.annotations_dir,
                                                 'annotations.tsv')

        # 9. GEN: <root>/antismash / antismash_dir=<override>
        self.antismash_dir = self._overrides.get(
            self.OR_ANTISMASH) or os.path.join(self._root, 'antismash')
        self.antismash_cache = {}

        # 10. GEN: <root>/bigscape / bigscape_dir=<override>
        self.bigscape_dir = self._overrides.get(
            self.OR_BIGSCAPE) or os.path.join(self._root, 'bigscape')
        # what we really want here is the subdirectory containing the network/annotation files,
        # but in case this is the path to the top level bigscape output, try to figure it out automatically
        if not os.path.exists(os.path.join(self.bigscape_dir, 'NRPS')):
            new_bigscape_dir = find_bigscape_dir(self.bigscape_dir)
            if new_bigscape_dir is not None:
                logger.info(
                    'Updating bigscape_dir to discovered location {}'.format(
                        new_bigscape_dir))
                self.bigscape_dir = new_bigscape_dir

        # 11. GEN: <root>/mibig_json / mibig_json_dir=<override>
        self.mibig_json_dir = self._overrides.get(
            self.OR_MIBIG_JSON) or os.path.join(self._root, 'mibig_json')

        # 12. MISC: <root>/params.xml
        self.params_file = os.path.join(self._root, 'params.xml')

        # 13. MISC: <root>/description.txt
        self.description_file = os.path.join(self._root, 'description.txt')

        # 14. MISC: <root>/include_strains.csv / include_strains_file=<override>
        self.include_strains_file = self._overrides.get(
            self.OR_INCLUDE_STRAINS) or os.path.join(self._root,
                                                     'include_strains.csv')

        # 15. CLASS: <root>/canopus / canopus_dir=<override>
        self.canopus_dir = self._overrides.get(
            self.OR_CANOPUS) or os.path.join(self._root, 'canopus')

        # 15. CLASS: <root>/canopus / canopus_dir=<override>
        self.molnetenhancer_dir = self._overrides.get(
            self.OR_MOLNETENHANCER) or os.path.join(self._root,
                                                    'molnetenhancer')

        for f in self.required_paths():
            if not os.path.exists(f):
                raise FileNotFoundError(
                    'File/directory "{}" does not exist or is not readable!'.
                    format(f))

        for f in self.optional_paths():
            if not os.path.exists(f):
                logger.warning(
                    'Optional file/directory "{}" does not exist or is not readable!'
                    .format(f))

    def load(self, met_only):
        # load strain mappings first
        if not self._load_strain_mappings():
            return False

        if not self._load_metabolomics():
            return False

        if not met_only and not self._load_genomics():
            return False

        if not self._load_class_info():
            return False

        self._load_optional()

        # Restrict strain list to only relevant strains (those that are present
        # in both genomic and
        if not met_only:
            # TODO add a config file option for this?
            self._filter_only_common_strains()

            # if the user specified a set of strains to be explicitly included, filter
            # out everything except those strains
            self._filter_user_strains()

        # if we don't have at least *some* strains here it probably means missing mappings
        # or a complete failure to parse things, so bail out
        if len(self.strains) == 0:
            raise Exception(
                'Failed to find *ANY* strains, missing strain_mappings.csv?')

        return True

    def _filter_user_strains(self):
        """
        If the user has supplied a list of strains to be explicitly included, go through the
        existing sets of objects we have and remove any that only include other strains. This
        involves an initial round of removing BGC and Spectrum objects, then a further round
        of removing now-empty GCF and MolFam objects.
        """
        if len(self.include_only_strains) == 0:
            logger.info('No further strain filtering to apply')
            return

        logger.info(
            'Found a list of {} strains to retain, filtering objects'.format(
                len(self.include_only_strains)))

        # filter the main list of strains
        self.strains.filter(self.include_only_strains)

        if len(self.strains) == 0:
            logger.error(
                'Strain list has been filtered down until it is empty! ')
            logger.error(
                'This probably indicates that you tried to specifically include a set of strains that had no overlap with the set common to metabolomics and genomics data (see the common_strains.csv in the dataset folder for a list of these'
            )
            raise Exception(
                'No strains left after filtering, cannot continue!')

        # get the list of BGCs which have a strain found in the set we were given
        bgcs_to_retain = {
            bgc for bgc in self.bgcs if bgc.strain in self.include_only_strains
        }
        # get the list of spectra which have at least one strain in the set
        spectra_to_retain = {
            spec for spec in self.spectra for sstrain in spec.strains
            if sstrain in self.include_only_strains
        }

        logger.info('Current / filtered BGC counts: {} / {}'.format(
            len(self.bgcs), len(bgcs_to_retain)))
        logger.info('Current / filtered spectra counts: {} / {}'.format(
            len(self.spectra), len(spectra_to_retain)))

        self.bgcs = list(bgcs_to_retain)
        for i, bgc in enumerate(self.bgcs):
            bgc.id = i

        self.spectra = list(spectra_to_retain)
        # also need to filter the set of strains attached to each spectrum
        for i, spec in enumerate(self.spectra):
            spec.strains.filter(self.include_only_strains)
            spec.id = i

        # now filter GCFs and MolFams based on the filtered BGCs and Spectra
        gcfs = {parent for bgc in self.bgcs for parent in bgc.parents}
        logger.info('Current / filtered GCF counts: {} / {}'.format(
            len(self.gcfs), len(gcfs)))
        self.gcfs = list(gcfs)
        # filter each GCF's strain list
        for i, gcf in enumerate(self.gcfs):
            gcf.strains.filter(self.include_only_strains)
            gcf.id = i

        molfams = {spec.family for spec in self.spectra}
        logger.info('Current / filtered MolFam counts: {} / {}'.format(
            len(self.molfams), len(molfams)))
        self.molfams = list(molfams)
        for i, molfam in enumerate(self.molfams):
            molfam.id = i

    def _filter_only_common_strains(self):
        """
        Filter strain population to only strains present in both genomic and molecular data
        """
        # TODO: Maybe there should be an option to specify which strains are used, both so we can
        #    selectively exclude strains, and include strains that are missing from either side.
        bgc_strains = {x.strain for x in self.bgcs}
        spectrum_strains = set().union(*[x.strains for x in self.spectra])
        common_strains = bgc_strains.intersection(spectrum_strains)
        logger.debug(
            'Filtering strains: genomics count {}, metabolomics count: {}'.
            format(len(bgc_strains), len(spectrum_strains)))
        logger.debug(f'Common strains found: {len(common_strains)}')

        # write out a list of the common strains to the dataset folder (might be useful for
        # anyone wanting to do additional filtering)
        cs_path = os.path.join(self._root, 'common_strains.csv')
        logger.info(f'Writing common strain labels to {cs_path}')
        with open(cs_path, 'w') as cs:
            cs.write('# strain label\n')
            for strain in self.strains:
                cs.write(f'{strain.id}\n')

        # filter the master list of strains down to include only the common set
        self.strains.filter(common_strains)

        for gcf in self.gcfs:
            gcf.strains.filter(common_strains)
        for spec in self.spectra:
            spec.strains.filter(common_strains)
        logger.info('Strains filtered down to total of {}'.format(
            len(self.strains)))

    def _load_optional(self):
        self.gnps_params = {}
        if os.path.exists(self.params_file):
            logger.debug('Loading params.xml')
            tree = ET.parse(self.params_file)
            root = tree.getroot()
            # this file has a simple structure:
            # <parameters>
            #   <parameter name="something">value</parameter>
            # </parameters>
            for param in root:
                self.gnps_params[param.attrib['name']] = param.text

            logger.debug(f'Parsed {len(self.gnps_params)} GNPS params')

        self.description_text = '<no description>'
        if os.path.exists(self.description_file):
            self.description_text = open(self.description_file).read()
            logger.debug('Parsed description text')

        self.include_only_strains = set()
        if os.path.exists(self.include_strains_file):
            logger.debug('Loading include_strains from {}'.format(
                self.include_strains_file))
            strain_list = open(self.include_strains_file).readlines()
            self.include_only_strains = StrainCollection()
            for line_num, sid in enumerate(strain_list):
                sid = sid.strip()  # get rid of newline
                strain_obj = self.strains.lookup(sid)
                if strain_obj is None:
                    logger.warning(
                        'Line {} of {}: invalid/unknown strain ID "{}"'.format(
                            line_num + 1, self.include_strains_file, sid))
                    continue
                self.include_only_strains.add(strain_obj)
            logger.debug('Found {} strain IDs in include_strains'.format(
                len(self.include_only_strains)))

    def _load_genomics_extra(self):
        if not os.path.exists(self.mibig_json_dir):
            if self._use_mibig:
                logger.info(
                    'Attempting to download MiBIG JSON database (v{})...'.
                    format(self._mibig_version))
                download_and_extract_mibig_metadata(self._root,
                                                self.mibig_json_dir,
                                                self._mibig_version)
            else:
                logger.warning(
                    'Not downloading MiBIG database automatically, use_mibig = false'
                )

        if not os.path.exists(self.bigscape_dir):
            should_run_bigscape = self._docker.get('run_bigscape',
                                                   self.RUN_BIGSCAPE_DEFAULT)
            extra_bigscape_parameters = self._docker.get(
                'extra_bigscape_parameters',
                self.EXTRA_BIGSCAPE_PARAMS_DEFAULT)
            if should_run_bigscape:
                # TODO this should not be attempted if not in Docker env
                logger.info(
                    'Running BiG-SCAPE! extra_bigscape_parameters="{}"'.format(
                        extra_bigscape_parameters))
                try:
                    run_bigscape('bigscape.py',
                                 os.path.join(self._root, 'antismash'),
                                 os.path.join(self._root, 'bigscape'),
                                 self.PFAM_PATH,
                                 extra_params=extra_bigscape_parameters)
                except Exception as e:
                    logger.warning(
                        'Failed to run BiG-SCAPE on antismash data, error was "{}"'
                        .format(e))

                self.bigscape_dir = find_bigscape_dir(self.bigscape_dir)

    # CG: load gemonics data into memory
    def _load_genomics(self):
        # TODO split this method up a bit

        # hmmscan apparently can't cope with filenames that have spaces in them, and so if you try
        # to run bigscape on a dataset like that it will just refuse. so, to try and avoid anyone
        # having to manually fix this stupid problem, we will attempt to rename every .gbk file here
        # TODO is this something nplinker should do, or dataset authors??
        logger.debug("\nLoading genomics starts...")
        logger.debug('Collecting .gbk files (and possibly renaming)')
        gbk_files = []
        t = time.time()
        renamed = 0

        # if the option is enabled, check for spaces in the folder names under
        # <dataset>/antismash and rename them, replacing spaces with underscores
        ignore_spaces = self._antismash_ignore_spaces
        if not ignore_spaces:
            logger.debug('Checking for spaces in antiSMASH folder names...')
            for root, dirs, files in os.walk(self.antismash_dir):
                for d in dirs:
                    if d.find(' ') != -1:
                        os.rename(os.path.join(root, d),
                                  os.path.join(root, d.replace(' ', '_')))
                        logger.warn(
                            'Renaming antiSMASH folder "{}" to "{}" to remove spaces! (suppress with ignore_spaces = true in config file)'
                            .format(d, d.replace(' ', '_')))

        # now want to gather a list of all .gbk files, and if necessary+enabled
        # rename them to replace spaces with underscores
        for root, dirs, files in os.walk(self.antismash_dir):
            for f in files:
                if f.lower().endswith('.gbk'):
                    fullpath = os.path.join(root, f)

                    if not ignore_spaces and fullpath.find(' ') != -1:
                        logger.debug(
                            'Spaces found in antiSMASH file "{}", renaming it'.
                            format(fullpath))
                        newpath = fullpath.replace(' ', '_')
                        os.rename(fullpath, newpath)
                        fullpath = newpath
                        renamed += 1

                    gbk_files.append(fullpath)

        logger.debug(f'.gbk collection took {time.time() - t:.3f}s')
        if renamed > 0:
            logger.info(
                '{}/{} .gbk files were renamed to remove spaces!'.format(
                    renamed, len(gbk_files)))

        # both the bigscape and mibig_json dirs expected by nplinker may not exist at this point. in some
        # cases this will cause an error later in the process, but this can also potentially be
        # resolved automatically:
        #   mibig_json => download and extract the JSON database
        #   bigscape => run BiG-SCAPE before continuing (if using the Docker image)
        self._load_genomics_extra()


        logger.debug(f'MibigBGCLoader({self.mibig_json_dir})')
        mibig_bgc_loader = MibigBGCLoader(self.mibig_json_dir)
        self.mibig_bgc_dict = mibig_bgc_loader.bgcs

        # add mibig bgc strains
        for bgc in self.mibig_bgc_dict.values():
            self.strains.add(bgc.strain)

        logger.debug('mibig_bgc_dict has {} entries'.format(
            len(self.mibig_bgc_dict)))

        missing_anno_files, missing_cluster_files, missing_network_files = [], [], []
        cluster_files, anno_files, network_files = {}, {}, {}

        # CG: bigscape data used in NPLinker
        # Take chemical class PKSI as example:
        #   Network_Annotations_PKSI.tsv
        #   PKSI_c0.30.network.tsv
        #   PKSI_clustering_c0.30.tsv

        for folder in self.BIGSCAPE_PRODUCT_TYPES:
            folder_path = os.path.join(self.bigscape_dir, folder)
            cutoff_filename = '{}_clustering_c0.{:02d}.tsv'.format(
                folder, self._bigscape_cutoff)
            cluster_filename = os.path.join(folder_path, cutoff_filename)
            annotation_filename = os.path.join(
                folder_path, f'Network_Annotations_{folder}.tsv')
            network_filename = os.path.join(
                folder_path,
                f'{folder}_c0.{self._bigscape_cutoff:02d}.network')

            # mandatory
            if not os.path.exists(annotation_filename):
                missing_anno_files.append(annotation_filename)
            else:
                anno_files[folder] = annotation_filename

            # also mandatory (?)
            if not os.path.exists(cluster_filename):
                missing_cluster_files.append(cluster_filename)
            else:
                cluster_files[folder] = cluster_filename

            # optional (?)
            if not os.path.exists(network_filename):
                missing_network_files.append(network_filename)
            else:
                network_files[folder] = network_filename

        # no files found here indicates a problem!
        if len(anno_files) == 0:
            raise Exception(
                'Failed to find *any* BiGSCAPE Network_Annotations tsv files under "{}" (incorrect cutoff value? currently set to {})'
                .format(self.bigscape_dir, self._bigscape_cutoff))

        def _list_missing_files(m, l):
            if len(l) == 0:
                return
            logger.warning(m)
            for i, fn in enumerate(l):
                logger.warning(f'  {i + 1}/{len(l)}: ' + fn)

        _list_missing_files(
            f'{len(missing_anno_files)} missing annotation tsv files:',
            missing_anno_files)
        _list_missing_files(
            '{} missing clustering tsv files:'.format(
                len(missing_cluster_files)), missing_cluster_files)
        _list_missing_files(
            f'{len(missing_network_files)} missing network files:',
            missing_network_files)

        # exclude any product types that don't have both annotation and cluster files
        self.product_types = []
        for prodtype in self.BIGSCAPE_PRODUCT_TYPES:
            if prodtype not in anno_files or prodtype not in cluster_files:
                # remove this product type completely
                for d in [anno_files, cluster_files, network_files]:
                    if prodtype in d:
                        del d[prodtype]
                logger.warning(
                    'Product type {} will be skipped due to missing files!'.
                    format(prodtype))
            else:
                self.product_types.append(prodtype)

        # generate a cache of antismash filenames to make matching them to BGC objects easier
        self.antismash_cache = {}
        logger.debug('Generating antiSMASH filename cache...')
        t = time.time()
        for f in gbk_files:
            path, filename = os.path.split(f)
            # cache with key set to bare filename with no extension
            basename = os.path.splitext(filename)[0]
            self.antismash_cache[basename] = f
            # also insert it with the folder name as matching on filename isn't always enough apparently
            parent = os.path.split(path)[-1]
            self.antismash_cache[f'{parent}_{basename}'] = f
        logger.debug(f'Cache generation took {time.time() - t:.3f}s')

        logger.debug(
            'loadBGC_from_cluster_files(antismash_dir={}, delimiters={})'.
            format(self.antismash_dir, self._antismash_delimiters))

        self.gcfs, self.bgcs, self.strains, unknown_strains = loadBGC_from_cluster_files(
            self.strains,
            cluster_files,
            anno_files,
            network_files,
            self.mibig_bgc_dict,
            self.mibig_json_dir,
            antismash_dir=self.antismash_dir,
            antismash_filenames=self.antismash_cache,
            antismash_format=self._antismash_format,
            antismash_delimiters=self._antismash_delimiters)

        us_path = os.path.join(self._root, 'unknown_strains_gen.csv')
        logger.warning(
            f'Writing unknown strains from GENOMICS data to {us_path}')
        with open(us_path, 'w') as us:
            us.write('# unknown strain label, filename\n')
            for strain, filename in unknown_strains.items():
                us.write(f'{strain}, {filename}\n')

        logger.debug("Loading genomics completed\n")

        return True

    def _load_metabolomics(self):
        spec_dict, self.spectra, self.molfams, unknown_strains = load_dataset(
            self.strains, self.mgf_file, self.edges_file, self.nodes_file,
            self.quantification_table_file, self.metadata_table_file,
            self._extended_metadata_table_parsing)

        us_path = os.path.join(self._root, 'unknown_strains_met.csv')
        logger.warning(
            'Writing unknown strains from METABOLOMICS data to {}'.format(
                us_path))
        with open(us_path, 'w') as us:
            us.write('# unknown strain label\n')
            for strain in unknown_strains.keys():
                us.write(f'{strain}\n')

        # load any available annotations from GNPS or user-provided files
        logger.info('Loading provided annotation files ({})'.format(
            self.annotations_dir))
        self.spectra = load_annotations(self.annotations_dir,
                                        self.annotations_config_file,
                                        self.spectra, spec_dict)
        return True

    def _load_class_info(self):
        """Load class match info (based on mibig) and chemical class predictions

        Run CANOPUS if asked for. First sirius is run through docker, if this
        fails, it is run with a version present on the path.

        Return:
            True if everything completes
        """
        # load Class_matches with mibig info from data
        mibig_class_file = self.datadir.joinpath(
            'MIBiG2.0_compounds_with_AS_BGC_CF_NPC_classes.txt')
        self.class_matches = ClassMatches(mibig_class_file)

        # run canopus if canopus_dir does not exist
        should_run_canopus = self._docker.get('run_canopus',
                                              self.RUN_CANOPUS_DEFAULT)
        extra_canopus_parameters = self._docker.get(
            'extra_canopus_parameters', self.EXTRA_CANOPUS_PARAMS_DEFAULT)
        if should_run_canopus:
            # don't run canopus when canopus dir exists already
            if not os.path.isdir(self.canopus_dir):
                logger.info(
                    'Running CANOPUS! extra_canopus_parameters="{}"'.format(
                        extra_canopus_parameters))
                try:
                    run_canopus(self.mgf_file, self.canopus_dir,
                                extra_canopus_parameters)
                except Exception as e:
                    logger.warning(
                        'Failed to run CANOPUS on mgf file with docker, error was "{}"'
                        .format(e))
                    logger.info(
                        'Trying to run CANOPUS again using SIRIUS from path')
                    try:
                        run_canopus(self.mgf_file, self.canopus_dir,
                                    extra_canopus_parameters)
                    except Exception as e:
                        logger.warning(
                            'Again failed to run CANOPUS on mgf file using sirius from path, error was "{}"'
                            .format(e))
            else:
                logger.info('Found CANOPUS dir, CANOPUS not run again!')

        # load Chem_class_predictions (canopus, molnetenhancer are loaded)
        chem_classes = ChemClassPredictions(self.canopus_dir,
                                            self.molnetenhancer_dir,
                                            self._root)
        # if no molfam classes transfer them from spectra (due to old style MN)
        if not chem_classes.canopus.molfam_classes and \
                chem_classes.canopus.spectra_classes:
            logger.debug("Added chemical compound classes for MFs")
            chem_classes.canopus.transfer_spec_classes_to_molfams(self.molfams)
        # include them in loader
        self.chem_classes = chem_classes
        return True

    def _load_strain_mappings(self):
        # this file should be a csv file, one line per strain, containing a list
        # of possible alternative IDs (the first one being the preferred ID).
        #
        # this is a per-dataset mapping, and is then merged with the global mapping file
        # packaged with nplinker itself
        self._init_global_strain_id_mapping()

        # now load the dataset mapping in the same way
        # TODO: what happens in case of clashes (differing primary IDs?)
        if not os.path.exists(self.strain_mappings_file):
            # create an empty placeholder file and show a warning
            logger.warn(
                'No strain_mappings.csv file found! Attempting to create one')
            self.strains.generate_strain_mappings(self.strain_mappings_file,
                                     self.antismash_dir)
            # raise Exception('Unable to load strain_mappings file: {}'.format(self.strain_mappings_file))
        else:
            self.strains.add_from_file(self.strain_mappings_file)
            logger.info('Loaded dataset strain IDs ({} total)'.format(
                len(self.strains)))

        return True

    def _init_global_strain_id_mapping(self):
        self.strains = StrainCollection()
        global_strain_id_file = self.datadir.joinpath('strain_id_mapping.csv')
        self.strains.add_from_file(global_strain_id_file)
        logger.info('Loaded global strain IDs ({} total)'.format(
            len(self.strains)))

    def required_paths(self):
        # these are files/paths that *must* exist for loading to begin
        return [
            self.nodes_file, self.edges_file, self.mgf_file, self.antismash_dir
        ]

    def optional_paths(self):
        return [self.annotations_dir]

    def webapp_scoring_cutoff(self):
        return self._webapp.get('tables_metcalf_threshold',
                                self.TABLES_CUTOFF_DEFAULT)

    def __repr__(self):
        return 'Root={}\n   MGF={}\n   EDGES={}\n   NODES={}\n   BIGSCAPE={}\n   ANTISMASH={}\n'.format(
            self._root, self.mgf_file, self.edges_file, self.nodes_file,
            self.bigscape_dir, self.antismash_dir)
