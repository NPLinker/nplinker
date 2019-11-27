import os
import glob
import time

from .metabolomics import load_spectra
from .metabolomics import load_edges
from .metabolomics import load_metadata
from .metabolomics import make_families

from .genomics import loadBGC_from_cluster_files
from .genomics import make_mibig_bgc_dict

from .annotations import load_annotations

from .strains import StrainCollection

from .logconfig import LogConfig
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
    OR_ANNO         = 'annotations_dir'
    OR_ANNO_CONFIG  = 'annotations_config_file'
    # and the same for genomics data
    OR_ANTISMASH    = 'antismash_dir'
    OR_BIGSCAPE     = 'bigscape_dir'
    OR_MIBIG_JSON   = 'mibig_json_dir'
    OR_STRAINS      = 'strain_mappings_file'

    BIGSCAPE_PRODUCT_TYPES = ['NRPS', 'Others', 'PKSI', 'PKS-NRP_Hybrids', 'PKSother', 'RiPPs', 'Saccharides', 'Terpene']

    def __init__(self, config_data):
        self._config = config_data
        self._dataset = config_data['dataset']
        self._overrides = self._dataset['overrides']
        self._antismash_format = self._dataset.get('antismash_format', self.ANTISMASH_FMT_DEFAULT)
        self._bigscape_cutoff = self._dataset.get('bigscape_cutoff', self.BIGSCAPE_CUTOFF_DEFAULT)
        self._root = self._config['dataset']['root']
        logger.debug('DatasetLoader({})'.format(self._root))

        # check antismash format is recognised
        if self._antismash_format not in self.ANTISMASH_FMTS:
            raise Exception('Unknown antismash format: {}'.format(self._antismash_format))

        # construct the paths and filenames required to load everything else and check 
        # they all seem to exist (but don't parse anything yet)

        # 1. strain mapping are used for everything else so
        self.strain_mappings_file = self._overrides.get(self.OR_STRAINS, os.path.join(self._root, 'strain_mappings.csv'))
        
        # 2. MET: <root>/clusterinfo_summary/<some UID>.tsv / nodes_file=<override>
        self.nodes_file = self._overrides.get(self.OR_NODES, find_via_glob(os.path.join(self._root, 'clusterinfo_summary', '*.tsv'), self.OR_NODES))

        # 3. MET: <root>/networkedges_selfloop/<some UID>.tsv / edges_file=<override>
        self.edges_file = self._overrides.get(self.OR_EDGES, find_via_glob(os.path.join(self._root, 'networkedges_selfloop', '*.selfloop'), self.OR_EDGES))

        # 4. MET: <root>/*.csv / extra_nodes_file=<override>
        # TODO is the glob input OK? 
        # => wait for updated dataset with latest output format
        self.extra_nodes_file = self._overrides.get(self.OR_EXTRA_NODES, find_via_glob(os.path.join(self._root, '*_quant.csv'), self.OR_EXTRA_NODES, optional=True))

        # 5. MET: <root>/spectra/specs_ms.mgf / mgf_file=<override>
        self.mgf_file = self._overrides.get(self.OR_MGF, os.path.join(self._root, 'spectra', 'specs_ms.mgf'))

        # 6. MET: <root>/metadata_table/metadata_table-<number>.txt / metadata_table_file=<override>
        self.metadata_table_file = self._overrides.get(self.OR_METADATA, find_via_glob(os.path.join(self._root, 'metadata_table', 'metadata_table*.txt'), self.OR_METADATA, optional=True))

        # 7. MET: <root>/quantification_table/quantification_table-<number>.csv / quantification_table_file=<override>
        self.quantification_table_file = self._overrides.get(self.OR_QUANT, find_via_glob(os.path.join(self._root, 'quantification_table', 'quantification_table*.csv'), self.OR_QUANT, optional=True))

        # 8. MET: <root>/DB_result/*.tsv / annotations_dir=<override>
        self.annotations_dir = self._overrides.get(self.OR_ANNO, os.path.join(self._root, 'DB_result'))
        self.annotations_config_file = self._overrides.get(self.OR_ANNO_CONFIG, os.path.join(self._root, 'DB_result', 'annotations.tsv'))
        
        # 9. GEN: <root>/antismash / antismash_dir=<override>
        self.antismash_dir = self._overrides.get(self.OR_ANTISMASH, os.path.join(self._root, 'antismash'))
        self.antismash_cache = {}

        # 10. GEN: <root>/bigscape / bigscape_dir=<override>
        self.bigscape_dir = self._overrides.get(self.OR_BIGSCAPE, os.path.join(self._root, 'bigscape'))

        # 11. GEN: <root>/mibig_json / mibig_json_dir=<override>
        self.mibig_json_dir = self._overrides.get(self.OR_MIBIG_JSON, os.path.join(self._root, 'mibig_json'))

        for f in self.required_paths():
            if not os.path.exists(f):
                raise FileNotFoundError('File/directory "{}" does not exist or is not readable!'.format(f))

        for f in self.optional_paths():
            if not os.path.exists(f):
                logger.warning('Optional file/directory "{}" does not exist or is not readable!'.format(f))

    def load(self):
        # load strain mappings first
        if not self._load_strain_mappings():
            return False

        if not self._load_metabolomics():
            return False

        if not self._load_genomics():
            return False

        return True

    def _load_genomics(self):
        logger.debug('make_mibig_bgc_dict({})'.format(self.mibig_json_dir))
        self.mibig_bgc_dict = make_mibig_bgc_dict(self.strains, self.mibig_json_dir)
        logger.debug('mibig_bgc_dict has {} entries'.format(len(self.mibig_bgc_dict)))

        missing_anno_files, missing_cluster_files, missing_network_files = [], [], []
        cluster_files, anno_files, network_files = {}, {}, {}

        for folder in self.BIGSCAPE_PRODUCT_TYPES:
            folder_path = os.path.join(self.bigscape_dir, folder)
            cutoff_filename = '{}_clustering_c0.{:02d}.tsv'.format(folder, self._bigscape_cutoff)
            cluster_filename = os.path.join(folder_path, cutoff_filename)
            annotation_filename = os.path.join(folder_path, 'Network_Annotations_{}.tsv'.format(folder))
            network_filename = os.path.join(folder_path, '{}_c0.{:02d}.network'.format(folder, self._bigscape_cutoff))

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
            raise Exception('Failed to find *any* BiGSCAPE Network_Annotations tsv files under "{}" (incorrect cutoff value? currently set to {})'.format(self.bigscape_dir, self._bigscape_cutoff))

        def _list_missing_files(t, l):
            if len(l) == 0:
                return
            logger.warning(t)
            for i, f in enumerate(l):
                logger.warning('  {}/{}: '.format(i+1, len(l)) + f)

        _list_missing_files('{} missing annotation tsv files:'.format(len(missing_anno_files)), missing_anno_files)
        _list_missing_files('{} missing clustering tsv files:'.format(len(missing_cluster_files)), missing_cluster_files)
        _list_missing_files('{} missing network files:'.format(len(missing_network_files)), missing_network_files)

        # exclude any product types that don't have both annotation and cluster files
        self.product_types = []
        for prodtype in self.BIGSCAPE_PRODUCT_TYPES:
            if prodtype not in anno_files or prodtype not in cluster_files:
                # remove this product type completely
                for d in [anno_files, cluster_files, network_files]:
                    if prodtype in d:
                        del d[prodtype]
                logger.warning('Product type {} will be skipped due to missing files!'.format(prodtype))
            else:
                self.product_types.append(prodtype)

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
                                                self.strains,
                                                cluster_files,
                                                anno_files,
                                                network_files,
                                                self.mibig_bgc_dict,
                                                antismash_dir=self.antismash_dir,
                                                antismash_filenames=self.antismash_cache,
                                                antismash_format=self._antismash_format)
        return True

    def _load_metabolomics(self):
        logger.debug('load_spectra({})'.format(self.mgf_file))
        self.spectra = load_spectra(self.mgf_file)

        spec_dict = {spec.spectrum_id : spec for spec in self.spectra}

        logger.debug('load_edges({})'.format(self.edges_file))
        self.spectra = load_edges(self.edges_file, self.spectra, spec_dict)

        logger.debug('load_metadata({})'.format(self.nodes_file))
        self.spectra = load_metadata(self.strains, self.nodes_file, self.extra_nodes_file, self.metadata_table_file, self.spectra, spec_dict)

        # load any available annotations from GNPS or user-provided files
        logger.info('Loading provided annotation files ({})'.format(self.annotations_dir))
        self.spectra = load_annotations(self.annotations_dir, self.annotations_config_file, self.spectra, spec_dict)

        self.molfams = make_families(self.spectra)
        logger.debug('make_families generated {} molfams'.format(len(self.molfams)))
        return True

    def _load_strain_mappings(self):
        # this file should be a csv file, one line per strain, containing a list
        # of possible alternative IDs (the first one being the preferred ID). 
        # 
        # this is a per-dataset mapping, and is then merged with the global mapping file
        # packaged with nplinker itself 
        self.strains = StrainCollection()

        global_strain_id_file = os.path.join(os.path.dirname(__file__), 'data', 'strain_id_mapping.csv')
        self.strains.add_from_file(global_strain_id_file)
        logger.info('Loaded global strain IDs ({} total)'.format(len(self.strains)))

        # now load the dataset mapping in the same way
        # TODO: what happens in case of clashes (differing primary IDs?)
        if not os.path.exists(self.strain_mappings_file):
            raise Exception('Unable to load strain_mappings file: {}'.format(self.strain_mappings_file))
        else:
            self.strains.add_from_file(self.strain_mappings_file)
            logger.info('Loaded dataset strain IDs ({} total)'.format(len(self.strains)))

        return True

    def required_paths(self):
        return [self.nodes_file, self.edges_file, self.mgf_file, self.antismash_dir, self.bigscape_dir, self.mibig_json_dir]

    def optional_paths(self):
        return [self.annotations_dir]

    def __repr__(self):
        return 'Root={}\n   MGF={}\n   EDGES={}\n   NODES={}\n   BIGSCAPE={}\n   ANTISMASH={}\n'.format(
                self._root, self.mgf_file, self.edges_file, self.nodes_file, self.bigscape_dir, self.antismash_dir)
