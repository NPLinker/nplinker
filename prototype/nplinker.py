import sys
import os
import glob
import logging

from metabolomics import load_spectra
from metabolomics import load_edges
from metabolomics import load_metadata
from metabolomics import make_families
from metabolomics import Spectrum
from metabolomics import MolecularFamily

from genomics import GCF
from genomics import loadBGC_from_cluster_files
from genomics import make_mibig_bgc_dict

from data_linking import DataLinks
from data_linking import RandomisedDataLinks
from data_linking import LinkLikelihood
from data_linking import LinkFinder

from logconfig import LogConfig

logger = LogConfig.getLogger(__file__)

class NPLinker(object):

    FOLDERS = ['NRPS', 'Others', 'PKSI', 'PKS-NRP_Hybrids', 'PKSother', 'RiPPs', 'Saccharides', 'Terpene']

    def __init__(self, mgf_file, edges_file, nodes_file, mibig_json_dir,
                    root_path, antismash_dir, antismash_format='flat'):

        for f in [mgf_file, edges_file, nodes_file, mibig_json_dir, root_path, antismash_dir]:
            if f is not None and not os.path.exists(f):
                raise FileNotFoundError('File/directory "{}" does not exist or is not readable!'.format(f))

        self.mgf_file = mgf_file
        self.edges_file = edges_file
        self.nodes_file = nodes_file
        self.mibig_json_dir = mibig_json_dir
        self.root_path = root_path
        self.antismash_dir = antismash_dir
        self.antismash_format = antismash_format

        self._spectra = []
        self._bgcs = []
        self._gcfs = []
        self._strains = []
        self._metadata = {}
        self._families = None
        self._mibig_bgc_dict = {}

        self._indices = {}
        self._edges = {}

        self._datalinks = None
        self._rdatalinks = []
        self._linkfinders = {}

    def load_data(self):
        logger.debug('load_data')

        logger.debug('load_spectra({})'.format(self.mgf_file))
        self._spectra = load_spectra(self.mgf_file)

        logger.debug('load_edges({})'.format(self.edges_file))
        load_edges(self._spectra, self.edges_file)

        # families = make_families(self.spectra) # TODO?
        logger.debug('load_metadata({})'.format(self.nodes_file))
        self._metadata = load_metadata(self._spectra, self.nodes_file)

        input_files, ann_files = [], []
        if self.mibig_json_dir is not None:
            logger.debug('make_mibig_bgc_dict({})'.format(self.mibig_json_dir))
            self._mibig_bgc_dict = make_mibig_bgc_dict(self.mibig_json_dir)

        for folder in NPLinker.FOLDERS:
            fam_file = os.path.join(self.root_path, folder)
            cluster_file = glob.glob(fam_file + os.sep + folder + "_clustering*")
            annotation_files = glob.glob(fam_file + os.sep + "Network_*")
            input_files.append(cluster_file[0])
            ann_files.append(annotation_files[0])

        logger.debug('loadBGC_from_cluster_files(antismash_dir={})'.format(self.antismash_dir))
        self._gcfs, self._bgcs, self._strains = loadBGC_from_cluster_files(
                                                input_files,
                                                ann_files,
                                                antismash_dir=self.antismash_dir, antismash_format=self.antismash_format,
                                                mibig_bgc_dict=self._mibig_bgc_dict)

        self._indices = {
            GCF : {gcf: x for x, gcf in enumerate(self._gcfs)},
            Spectrum : {spectra: x for x, spectra in enumerate(self._spectra)},
            str : {strain: x for x, strain in enumerate(self._strains)},
        }

        logger.debug('load_data: completed')
        return True

    def process_dataset(self, find_correlations=True, random_count=50,
                            scoring_methods=['metcalf', 'likescore']):
        """
        Construct the DataLinks and LinkFinder objects from the loaded dataset.
        Then run scoring functions using the methods in scoring_methods with
        every created DataLinks instance.
        """
        if len(self._spectra) == 0 or len(self._gcfs) == 0 or len(self._strains) == 0:
            logger.info('process_dataset: calling load_data')
            self.load_data()

        logger.debug('Creating a DataLinks object')
        self._datalinks = DataLinks()
        self._datalinks.load_data(self._spectra, self._gcfs, self._strains)
        logger.debug('DataLinks load_data complete')
        if find_correlations:
            logger.debug('Finding correlations')
            self._datalinks.find_correlations()

        logger.debug('Generating RandomisedDataLinks, n={}'.format(random_count))
        self._rdatalinks = [RandomisedDataLinks.from_datalinks(self._datalinks, find_correlations) for x in range(random_count)]

        # TODO would it make more sense to make each LinkFinder encapsulate a DataLinks object?
        self._linkfinders[self._datalinks] = LinkFinder()
        for rd in self._rdatalinks:
            self._linkfinders[rd] = LinkFinder()

        logger.info('Generating scores, enabled methods={}'.format(scoring_methods))

        for dl, lf in self._linkfinders.items():
            if 'metcalf' in scoring_methods:
                lf.metcalf_scoring(dl, type='spec-gcf')
                lf.metcalf_scoring(dl, type='fam-gcf')

            if 'likescore' in scoring_methods:
                ll = LinkLikelihood()
                ll.calculate_likelihoods(dl, type='spec-gcf')
                ll.calculate_likelihoods(dl, type='fam-gcf')
                lf.likelihood_scoring(dl, ll, type='spec-gcf')
                lf.likelihood_scoring(dl, ll, type='fam-gcf')
        
        return True

    def _get_ids_for_objects(self, objects):
        if isinstance(objects, list) and type(objects[0]) in self._indices:
            t = type(objects[0])
            return [self._indices[t][x] for x in objects]
        
        if type(objects) in self._indices:
            t = type(objects)
            return self._indices[t][objects]

        logger.warn('Unknown type supplied to _get_ids_for_objects ({})'.format(type(objects)))
        return None

    def get_links(self, objects, datalinks=None, scoring_method='metcalf', score_cutoff=0.5):
        if len(self._linkfinders) == 0:
            raise Exception('Need to call process_dataset first')

        obj_classes = [Spectrum, MolecularFamily, GCF]
        obj_types = ['spec', 'fam', 'gcf']
        obj_type = None

        object = objects[0] if isinstance(objects, list) else objects

        for t in range(len(obj_types)):
            if isinstance(object, obj_classes[t]):
                obj_type = obj_types[t]
                break

        if obj_type is None:
            given_type = type(objects[0]) if isinstance(objects, list) else type(objects)
            raise Exception('Unsupported type passed to get_scores ({})'.format(given_type))

        datalinks = datalinks or self._datalinks

        # TODO this will need updated when current implementation is changed to 
        # avoid using IDs
        objects = self._get_ids_for_objects(objects)
        if objects is None:
            raise Exception('Error mapping objects to IDs!')

        result = self._linkfinders[datalinks].get_links(objects, obj_type, scoring_method, score_cutoff)
        print(result, type(result))
        if isinstance(result, list):
            print(len(result))
            print(result[0], result[0].shape)
        return result

    @property
    def strains(self):
        return self._strains

    @property
    def gcfs(self):
        return self._gcfs

    @property
    def spectra(self):
        return self._spectra

    @property
    def metadata(self):
        return self._metadata

    @property
    def mibig_bgc_dict(self):
        return self._mibig_bgc_dict

    @property
    def datalinks(self):
        return self._datalinks

    @property
    def rdatalinks(self):
        return self._rdatalinks

if __name__ == "__main__":
    # TODO will need to handle arguments here
    # TODO possible to auto-discover some/all of these filenames for a given dataset?

    # crusemann dataset
    DATASET = '/mnt/archive/nplinker_data/crusemann'
    MGF_FILE = os.path.join(DATASET, 'gnps/METABOLOMICS-SNETS-c36f90ba-download_clustered_spectra-main.mgf')
    NODES_FILE = os.path.join(DATASET, 'gnps/0d51c5b6c73b489185a5503d319977ab..out')
    EDGES_FILE = os.path.join(DATASET, 'gnps/9a93d720f69143bb9f971db39b5d2ba2.pairsinfo')
    ROOT_PATH = os.path.join(DATASET, 'bigscape/bigscape_corason_crusemann_complete_annotated_mibigs_mix_automode_20180713/network_files/2018-07-13_16-34-11_hybrids_auto_crusemann_bgcs_automode_mix_mibig')
    ANTISMASH_DIR = os.path.join(DATASET, 'antismash/justin-20181022/')

    # TODO is this optional???
    # MIBIG_JSON_DIR = DATASET + "\\Data\\mibig\\mibig_json_1.4"
    MIBIG_JSON_DIR = None

    # can set default logging configuration this way...
    LogConfig.setLogLevel(logging.DEBUG)

    # create top level "app" object passing in dataset details 
    # TODO alternative ways to specify this, e.g. in a .json file etc
    npl = NPLinker(MGF_FILE, EDGES_FILE, NODES_FILE, MIBIG_JSON_DIR, ROOT_PATH, ANTISMASH_DIR)

    # load the dataset
    if not npl.load_data():
        print('Failed to load the dataset!')
        sys.exit(-1)

    # generate all datalinks objects (including list of random ones)
    if not npl.process_dataset(random_count=1):
        print('Failed to process dataset')
        sys.exit(-1)

    # can now call get_links to get edge info
    gcfs = npl.gcfs[:15]
    gcf_links = npl.get_links(gcfs)

    # or do the same with a custom datalinks object
    r_gcf_links = npl.get_links(gcfs, datalinks=npl.rdatalinks[0])

    # or use a specific scoring method (defaults to metcalf)
    spectra = npl.spectra[100:115]
    spectra_links = npl.get_links(spectra, scoring_method='likescore')
    
