import os
import glob
import time

from metabolomics import load_spectra
from metabolomics import load_edges
from metabolomics import load_metadata
from metabolomics import make_families

from genomics import loadBGC_from_cluster_files
from genomics import make_mibig_bgc_dict

from logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

def find_via_glob(path, file_type, optional=False):
    try:
        filename = glob.glob(path)[0]
        return filename
    except (OSError, IndexError) as e:
        if not optional:
            # "from None" suppresses the traceback for the original exception, which isn't really needed
            raise Exception('ERROR: unable to find {} in path "{}"'.format(file_type, path)) from None

        logger.warn('WARNING: unable to find {} in path "{}"'.format(file_type, path))
        return None

class DatasetLoader(object):

    ANTISMASH_FMT_DEFAULT   = 'default'
    ANTISMASH_FMT_FLAT      = 'flat'
    ANTISMASH_FMTS          = [ANTISMASH_FMT_DEFAULT, ANTISMASH_FMT_FLAT]

    BIGSCAPE_CUTOFF_DEFAULT = 30

    # keys for overriding metabolomics data elements
    OR_NODES        = 'nodes_file'
    OR_EDGES        = 'edges_file'
    OR_EXTRA_NODES  = 'extra_nodes_file'
    OR_MGF          = 'mgf_file'
    OR_METADATA     = 'metadata_table_file'
    OR_QUANT        = 'quantification_table_file'
    # and the same for genomics data
    OR_ANTISMASH    = 'antismash_dir'
    OR_BIGSCAPE     = 'bigscape_dir'
    OR_MIBIG_JSON   = 'mibig_json_dir'

    # bigscape classes
    BIGSCAPE_CLASSES = ['NRPS', 'Others', 'PKSI', 'PKS-NRP_Hybrids', 'PKSother', 'RiPPs', 'Saccharides', 'Terpene']

    def __init__(self, config_data):
        self._config = config_data
        self._dataset = config_data['dataset']
        self._overrides = self._dataset['overrides']
        self._antismash_format = self._dataset.get('antismash_format', self.ANTISMASH_FMT_DEFAULT)
        self._bigscape_cutoff = self._dataset.get('bigscape_cutoff', self.BIGSCAPE_CUTOFF_DEFAULT)
        self._root = self._config['dataset']['root']
        logger.debug('DatasetLoader({})'.format(self._root))

        # check antismash format is recognised

        # construct the paths and filenames required to load everything and check 
        # they all seem to exist (but don't parse anything yet)
        
        # 1. MET: <root>/clusterinfo_summary/<some UID>.tsv / nodes_file=<override>
        self.nodes_file = self._overrides.get(self.OR_NODES, find_via_glob(os.path.join(self._root, 'clusterinfo_summary', '*.tsv'), self.OR_NODES))

        # 2. MET: <root>/networkedges_selfloop/<some UID>.tsv / edges_file=<override>
        self.edges_file = self._overrides.get(self.OR_EDGES, find_via_glob(os.path.join(self._root, 'networkedges_selfloop', '*.selfloop'), self.OR_EDGES))

        # 3. MET: <root>/*.csv / extra_nodes_file=<override>
        self.extra_nodes_file = self._overrides.get(self.OR_EXTRA_NODES, find_via_glob(os.path.join(self._root, '*.csv'), self.OR_EXTRA_NODES, optional=True))

        # 4. MET: <root>/spectra/specs_ms.mgf / mgf_file=<override>
        self.mgf_file = self._overrides.get(self.OR_MGF, os.path.join(self._root, 'spectra', 'specs_ms.mgf'))

        # 5. MET: <root>/metadata_table/metadata_table-<number>.txt / metadata_table_file=<override>
        self.metadata_table_file = self._overrides.get(self.OR_METADATA, find_via_glob(os.path.join(self._root, 'metadata_table', 'metadata_table*.txt'), self.OR_METADATA, optional=True))

        # 6. MET: <root>/quantification_table/quantification_table-<number>.csv / quantification_table_file=<override>
        self.quantification_table_file = self._overrides.get(self.OR_QUANT, find_via_glob(os.path.join(self._root, 'quantification_table', 'quantification_table*.csv'), self.OR_QUANT, optional=True))

        # 7. GEN: <root>/antismash / antismash_dir=<override>
        self.antismash_dir = self._overrides.get(self.OR_ANTISMASH, os.path.join(self._root, 'antismash'))
        self.antismash_cache = {}

        # 8. GEN: <root>/bigscape / bigscape_dir=<override>
        self.bigscape_dir = self._overrides.get(self.OR_BIGSCAPE, os.path.join(self._root, 'bigscape'))

        # 9. GEN: <root>/mibig_json / mibig_json_dir=<override>
        self.mibig_json_dir = self._overrides.get(self.OR_MIBIG_JSON, os.path.join(self._root, 'mibig_json'))

        for f in self.required_paths():
            if not os.path.exists(f):
                raise FileNotFoundError('File/directory "{}" does not exist or is not readable!'.format(f))

        # this is optional so just warn 
        if not os.path.exists(self.mibig_json_dir):
            logger.warning('mibig_json_dir "{}" does not exist or is not readable!'.format(self.mibig_json_dir))

    def load(self):
        # metabolomics stuff first
        logger.debug('load_spectra({})'.format(self.mgf_file))
        self.spectra = load_spectra(self.mgf_file)

        spec_dict = {spec.spectrum_id : spec for spec in self.spectra}

        logger.debug('load_edges({})'.format(self.edges_file))
        self.spectra = load_edges(self.edges_file, self.spectra, spec_dict)

        logger.debug('load_metadata({})'.format(self.nodes_file))
        self.spectra, self.strains = load_metadata(self.nodes_file, self.extra_nodes_file, self.spectra, spec_dict)

        self.molfams = make_families(self.spectra) # TODO?
        logger.debug('make_families generated {} molfams'.format(len(self.molfams)))

        # and now genomics
        cluster_files, anno_files = [], []
        logger.debug('make_mibig_bgc_dict({})'.format(self.mibig_json_dir))
        self.mibig_bgc_dict = make_mibig_bgc_dict(self.mibig_json_dir)
        logger.debug('mibig_bgc_dict has {} entries'.format(len(self.mibig_bgc_dict)))

        # TODO can this just be modified to search all folders in the root path instead of hardcoding them?
        for folder in self.BIGSCAPE_CLASSES:
            folder_path = os.path.join(self.bigscape_dir, folder)
            cutoff_filename = '{}_clustering_c0.{:02d}.tsv'.format(folder, self._bigscape_cutoff)
            cluster_filename = os.path.join(folder_path, cutoff_filename)
            annotation_filename = os.path.join(folder_path, 'Network_Annotations_{}.tsv'.format(folder))

            if not os.path.exists(annotation_filename) or not os.path.exists(cluster_filename):
                # TODO which folders are supposed to exist? is it a critical error if some of them
                # don't appear in a dataset???
                logger.warning('Missing annotation/cluster tsv files in folder: {}'.format(folder_path))
                continue

            cluster_files.append(cluster_filename)
            anno_files.append(annotation_filename)

        # no files found here indicates a problem!
        if len(anno_files) == 0:
            raise Exception('Failed to find any BIGSCAPE files under "{}" (incorrect cutoff value? currently set to {})'.format(self.bigscape_dir, self._bigscape_cutoff))

        # generate a cache of antismash filenames to make matching them to BGC objects easier
        self.antismash_cache = {}
        logger.debug('Generating antiSMASH filename cache...')
        t = time.time()
        for root, dirs, files in os.walk(self.antismash_dir):
            for f in files:
                if f.endswith('.gbk'):
                    basename = os.path.splitext(f)[0]
                    fullpath = os.path.join(root, f)
                    self.antismash_cache[basename] = fullpath
                    # also insert it with the folder name as matching on filename isn't always enough apparently
                    parent = os.path.split(root)[-1]
                    self.antismash_cache['{}_{}'.format(parent, basename)] = fullpath
        logger.debug('Cache generation took {:.3f}s'.format(time.time() - t))

        logger.debug('loadBGC_from_cluster_files(antismash_dir={})'.format(self.antismash_dir))
        self.gcfs, self.bgcs, self.strains = loadBGC_from_cluster_files(
                                                cluster_files,
                                                anno_files,
                                                antismash_dir=self.antismash_dir,
                                                antismash_filenames=self.antismash_cache,
                                                antismash_format=self._antismash_format,
                                                mibig_bgc_dict=self.mibig_bgc_dict)

        return True

    def required_paths(self):
        # return [self.nodes_file, self.edges_file, self.extra_nodes_file, self.mgf_file, self.metadata_table_file, self.quantification_table_file, self.antismash_dir, self.bigscape_dir]
        return [self.nodes_file, self.edges_file, self.mgf_file, self.antismash_dir, self.bigscape_dir]

    def __repr__(self):
        return 'Root={}\n   MGF={}\n   EDGES={}\n   NODES={}\n   BIGSCAPE={}\n   ANTISMASH={}\n'.format(
                self._root, self.mgf_file, self.edges_file, self.nodes_file, self.bigscape_dir, self.antismash_dir)
