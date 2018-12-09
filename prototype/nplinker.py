import sys
import os
import glob
import logging
import pickle

import numpy as np

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
    SCORING_METHODS = ['metcalf', 'likescore']
    DEFAULT_SCORING = 'metcalf'
    OBJ_CLASSES = [Spectrum, MolecularFamily, GCF]

    def __init__(self, mgf_file, edges_file, nodes_file, mibig_json_dir,
                    root_path, antismash_dir, antismash_format='flat',
                    repro_file=None):

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

        self._spec_gcf_scores = {m: None for m in NPLinker.SCORING_METHODS}
        self._fam_gcf_scores = {m: None for m in NPLinker.SCORING_METHODS}

        self._datalinks = None
        self._rdatalinks = []
        self._linkfinders = {}

        self._repro_data = {}
        if repro_file is not None:
            self._save_repro_data(repro_file)

    def _collect_repro_data(self):
        """
        This method creates a dict containing various bits of information about
        the current execution of nplinker. This data can optionally be saved to
        a file in order to aid reproducibility.
        """

        self._repro_data = {}
        # TODO best way to embed latest git commit hash? probably with a packaging script...
        # TODO versions of all Python dependencies used (just rely on Pipfile.lock here?)
        
        # insert command line arguments
        self._repro_data['args'] = {}
        for i, arg in enumerate(sys.argv):
            self._repro_data['args'][i] = arg

        # TODO anything else?

    def _save_repro_data(self, filename):
        self._collect_repro_data()
        with open(filename, 'wb') as repro_file:
            # TODO better format than pickle? 
            pickle.dump(self._repro_data, repro_file)
            logger.info('Saving reproducibility data to {}'.format(filename))

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

        logger.debug('load_data: completed')
        return True

    def process_dataset(self, find_correlations=True, random_count=50,
                            scoring_methods=SCORING_METHODS):
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

        # for storing results of calls to get_links
        self._spec_gcf_scores = {m: np.zeros((len(self._spectra), len(self._gcfs))) for m in NPLinker.SCORING_METHODS}
        self._fam_gcf_scores = {m: np.zeros((self._datalinks.M_fam_strain.shape[0], len(self._gcfs))) for m in NPLinker.SCORING_METHODS}

        return True

    def get_links(self, objects, datalinks=None, scoring_method=DEFAULT_SCORING, scoring_cutoff=0.5):
        """
        This method is used to collect scores for objects of type GCF/Spectrum/MolecularFamily given
        a scoring method and cutoff value.

        The objects parameter can be a single instance of any of those 3 types, a list of instances
        of any one type, or a mixed list containing instances of all 3 types. In the first two cases,
        the return value will be a list containing instances of objects for which links/edges were
        found with scores above the supplied cutoff value. For the latter case, the return value
        will be a 3-element list containing the results for each of the 3 object types, always in
        the order [GCF, Spectrum, MolecularFamily]. 

        To look up the actual scores for the returned objects, use NPLinker.<method>_scores
        """
        # single object, pass it straight through
        if not isinstance(objects, list):
            logger.debug('get_links: single object')
            return self._get_links([objects], datalinks, scoring_method, scoring_cutoff)

        # check if the list contains one or multiple types of object
        type_counts = {NPLinker.OBJ_CLASSES[t]: 0 for t in range(len(NPLinker.OBJ_CLASSES))}
        for o in objects:
            try:
                type_counts[type(o)] += 1
            except KeyError:
                raise Exception('Unsupported type passed to get_links ({})'.format(type(o)))

        if max(type_counts.values()) == len(objects):
            # if the list seems to have a uniform type, again pass straight through
            logger.debug('get_links: uniformly typed list ({})'.format(type(objects[0])))
            return self._get_links(objects, datalinks, scoring_method, scoring_cutoff)

        # finally, if it's a mix of objects, extract them into separate lists and call 
        # the internal method multiple times instead
        obj_lists = {t: [] for t in NPLinker.OBJ_CLASSES}
        for o in objects:
            obj_lists[type(o)].append(o)

        # keep consistent ordering here: GCF, Spectrum, MolecularFamily
        obj_lists = [obj_lists[GCF], obj_lists[Spectrum], obj_lists[MolecularFamily]]

        results = [[] for t in NPLinker.OBJ_CLASSES]
        logger.debug('get_links: multiple lists')
        for i in range(len(NPLinker.OBJ_CLASSES)):
            if len(obj_lists[i]) > 0:
                results[i].append(self._get_links(obj_lists[i], datalinks, scoring_method, scoring_cutoff))

        return results

    def _get_links(self, objects, datalinks, scoring_method, scoring_cutoff):
        if len(objects) == 0:
            return []

        if len(self._linkfinders) == 0:
            raise Exception('Need to call process_dataset first')

        if scoring_method not in NPLinker.SCORING_METHODS:
            raise Exception('unknown scoring method "{}"'.format(scoring_method))

        datalinks = datalinks or self._datalinks

        results = self._linkfinders[datalinks].get_links(datalinks, objects, scoring_method, scoring_cutoff)
        
        # first dimension of results contains input object ID, other object ID,
        # and the score for the link between them both
        R_SRC_ID, R_DST_ID, R_SCORE = range(3)
        scores_found = set()
        input_type = type(objects[0])

        if input_type == GCF:
            logger.debug('get_links: result type is GCF, inputs={}, results={}'.format(len(objects), len(results)))
            # for GCF input, results contains two arrays, giving spec-gcf and fam-gcf links respectively
            result_gcf_spec, result_gcf_fam = results[0], results[1]

            spec_scores = self._spec_gcf_scores[scoring_method]
            fam_scores = self._fam_gcf_scores[scoring_method]

            for res, sco, type_ in [(result_gcf_spec, spec_scores, 'spec'), (result_gcf_fam, fam_scores, 'fam')]:
                if res.shape[1] == 0:
                    logger.debug('Found no links for {}, [{}]'.format(objects, type_))
                    continue # no results

                logger.debug('Found {} links for {}, [{}]'.format(res.shape[1], objects, type_))

                # for each GCF in the result array
                for j, gcf_id in enumerate(res[R_SRC_ID]):
                    obj_id = int(res[R_DST_ID, j])
                    gcf_id = int(gcf_id)
                    # record the link between the Spectrum/MolecularFamily and this GCF
                    sco[obj_id, gcf_id] = res[R_SCORE, j]
                    logger.debug('Setting score ({}, {}) = {} [{}]'.format(obj_id, gcf_id, res[R_SCORE, j], type_))
                    scores_found.add(self._gcfs[gcf_id])
        else:
            logger.debug('get_links: result type is Spectrum/MolFam, inputs={}, results={}'.format(len(objects), len(results)))
            # for non-GCF input, result is a list containing a single array, shape (3, x)
            # where x is the number of links 
            results = results[0]
            if results.shape[1] == 0:
                logger.debug('Found no links for {}'.format(objects))
                return []

            obj_results = self._spec_gcf_scores[scoring_method] if input_type == Spectrum else self._fam_gcf_scores[scoring_method]
            logger.debug('Found {} links for {}'.format(results.shape[1], objects))

            # for each Spectrum/MolecularFamily in this result
            for i, obj_id in enumerate(results[R_SRC_ID]):
                gcf_id = int(results[R_DST_ID, i])
                obj_id = int(obj_id)
                # record the link between the GCF and this Spectrum/MolecularFamily
                obj_results[obj_id, gcf_id] = results[R_SCORE, i]
                logger.debug('Setting score ({}, {}) = {}'.format(obj_id, gcf_id, results[R_SCORE, i]))
                scores_found.add(self._spectra[obj_id] if input_type == Spectrum else self._families[obj_id])

        return list(scores_found)

    def get_common_strains(self, objects_a, objects_b):
        return self._datalinks.common_strains(objects_a, objects_b)

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
        """
        Returns the main DataLinks object used for scoring etc
        """
        return self._datalinks

    @property
    def rdatalinks(self):
        """
        Returns the list of RandomisedDataLinks objects
        """
        return self._rdatalinks

    @property
    def repro_data(self):
        """
        Returns the dict of reproducibility data
        """
        return self._repro_data

    # TODO better way of accessing the scores? or is this OK?
    @property
    def spec_scores(self, scoring_method=DEFAULT_SCORING):
        return self._spec_gcf_scores[scoring_method]
        
    @property
    def fam_scores(self, scoring_method=DEFAULT_SCORING):
        return self._fam_gcf_scores[scoring_method]

    @property
    def gcf_scores(self, scoring_method=DEFAULT_SCORING):
        return self._spec_gcf_scores[scoring_method], self._fam_gcf_scores[scoring_method]


if __name__ == "__main__":
    # TODO will need to handle arguments here (via argparse)
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

    # get_links with a single Spectrum/MolecularFamily
    obj = npl.spectra[2763]
    result = npl.get_links(obj, scoring_method='metcalf', scoring_cutoff=150)
    # if the input Spectrum had any links above the cutoff, it will appear in result:
    if len(result) == 1:
        print('Found links to object {}'.format(obj))
        # lookup the GCFs it has links to
        spec_scores = npl.spec_scores
        linked_gcfs = [npl.gcfs[gcf] for gcf in np.where(spec_scores[obj.id, :] > 150)[0]]

    # get_links with a list of GCFs
    objects = npl.gcfs[763:766]
    result = npl.get_links(objects, scoring_method='metcalf', scoring_cutoff=150)
    if len(result) > 0:
        print('{}/{} objects with links found!'.format(len(result), len(objects)))
        spec_scores = npl.spec_scores
        linked_spectra = {}
        for gcf in result:
            linked_spectra[gcf] = [npl.spectra[spec] for spec in np.where(spec_scores[:, gcf.id])[0]]
            print(gcf, linked_spectra[gcf], len(linked_spectra[gcf]))

    # finding common_strains
    common_strains = npl.get_common_strains(npl.gcfs[763:766], npl.spectra[3403])
    for spec, gcfs in common_strains.items():
        for gcf, strains in gcfs.items():
            for strain in strains:
                print('{} and {} share strain #{}'.format(spec, gcf, strain))

