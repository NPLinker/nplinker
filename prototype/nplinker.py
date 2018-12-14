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

# TODO Google-style docstrings

class NPLinker(object):

    FOLDERS = ['NRPS', 'Others', 'PKSI', 'PKS-NRP_Hybrids', 'PKSother', 'RiPPs', 'Saccharides', 'Terpene']
    SCORING_METHODS = ['metcalf', 'likescore']
    DEFAULT_SCORING = 'metcalf'
    OBJ_CLASSES = [Spectrum, MolecularFamily, GCF]
    # enumeration for accessing results of LinkFinder.get_links, which are (3, num_links) arrays:
    # - R_SRC_ID: the ID of an object that was supplied as input to get_links
    # - R_DST_ID: the ID of an object that was discovered to have a link to an input object
    # - R_SCORE: the score for the link between a pair of objects
    R_SRC_ID, R_DST_ID, R_SCORE = range(3)

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
        self._families = []
        self._mibig_bgc_dict = {}

        self._datalinks = None
        self._linkfinder = None

        self._random_scores = {}
        self._links = {}

        self._repro_data = {}
        if repro_file is not None:
            self._save_repro_data(repro_file)

    def _collect_repro_data(self):
        """Creates a dict containing data to aid reproducible runs of nplinker.

        This method creates a dict containing various bits of information about
        the current execution of nplinker. This data will typically be saved to
        a file in order to aid reproducibility using :func:`_save_repro_data`. 

        TODO describe contents

        Returns:
            A dict containing the information described above
        """

        self._repro_data = {}
        # TODO best way to embed latest git commit hash? probably with a packaging script...
        # TODO versions of all Python dependencies used (just rely on Pipfile.lock here?)
        
        # insert command line arguments
        self._repro_data['args'] = {}
        for i, arg in enumerate(sys.argv):
            self._repro_data['args'][i] = arg

        # TODO anything else to include here?

        return self._repro_data

    def _save_repro_data(self, filename):
        self._collect_repro_data()
        with open(filename, 'wb') as repro_file:
            # TODO is pickle the best format to use?
            pickle.dump(self._repro_data, repro_file)
            logger.info('Saving reproducibility data to {}'.format(filename))

    def load_data(self):
        """Loads the basic components of a dataset.

        This method is responsible for loading the various pieces of the supplied dataset into
        memory and doing any initial parsing/object manipulation required. After it completes,
        applications can access the lists of GCFs, Spectra, MolecularFamilies and strains 
        using the corresponding properties of the NPLinker class.

        Returns:
            True if successful, False otherwise
        """

        # TODO error handling (FileNotFoundError etc)
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

    def _generate_scores(self, dl, lf, scoring_methods):
        if 'metcalf' in scoring_methods:
            lf.metcalf_scoring(dl, type='spec-gcf')
            lf.metcalf_scoring(dl, type='fam-gcf')

        if 'likescore' in scoring_methods:
            ll = LinkLikelihood()
            ll.calculate_likelihoods(dl, type='spec-gcf')
            ll.calculate_likelihoods(dl, type='fam-gcf')
            lf.likelihood_scoring(dl, ll, type='spec-gcf')
            lf.likelihood_scoring(dl, ll, type='fam-gcf')

    def process_dataset(self, random_count=50, scoring_methods=SCORING_METHODS):
        """Construct the DataLinks and LinkFinder objects from the loaded dataset.

        Deals with initialising all the objects used for scoring/linking, and also
        currently manages the creation of randomised scoring matrices for the different
        object types. 
        """
        if len(self._spectra) == 0 or len(self._gcfs) == 0 or len(self._strains) == 0:
            logger.info('process_dataset: calling load_data')
            self.load_data()

        self._datalinks = DataLinks()
        self._datalinks.load_data(self._spectra, self._gcfs, self._strains)
        self._datalinks.find_correlations()
        self._linkfinder = LinkFinder()

        # generate the actual scores from the standard DataLinks object
        logger.info('Generating scores, enabled methods={}'.format(scoring_methods))
        self._generate_scores(self._datalinks, self._linkfinder, scoring_methods)

        # create the randomised scoring matrices, and stack them together for later use
        rdatalinks = [RandomisedDataLinks.from_datalinks(self._datalinks, True) for x in range(random_count)]
        rlinkfinders = [LinkFinder() for x in range(random_count)]
        logger.info('Generating randomised scores, enabled methods={}'.format(scoring_methods))
        
        for i in range(random_count):
            self._generate_scores(rdatalinks[i], rlinkfinders[i], scoring_methods)
        
        # create a set of dicts indexed by class to store the matrices
        self._random_scores = {}
        for m in scoring_methods:
            rand_scores = \
            {
                Spectrum:           np.zeros((len(self._spectra), 0)),
                MolecularFamily:    np.zeros((len(self._families), 0)),
                GCF:                {
                                        Spectrum:           np.zeros((0, len(self._gcfs))),
                                        MolecularFamily:    np.zeros((0, len(self._gcfs))),
                                    }
            }

            for i in range(random_count):
                rand_scores[GCF][Spectrum] = np.r_[rand_scores[GCF][Spectrum], rlinkfinders[i].get_scores(m, 'spec-gcf')]
                rand_scores[GCF][MolecularFamily] = np.r_[rand_scores[GCF][MolecularFamily], rlinkfinders[i].get_scores(m, 'fam-gcf')]
                rand_scores[Spectrum] = np.c_[rand_scores[Spectrum], rlinkfinders[i].get_scores(m, 'spec-gcf')]
                if len(self._families) > 0:
                    rand_scores[MolecularFamily] = np.c_[rand_scores[MolecularFamily], rlinkfinders[i].get_scores(m, 'fam-gcf')]

            self._random_scores[m] = rand_scores

        self.clear_links()
        return True

    def clear_links(self):
        """Clear any previously retrieved links"""
        self._links = {}

    def get_links(self, objects, datalinks=None, scoring_method=DEFAULT_SCORING, sig_percentile=95, scoring_cutoff=None):
        """Collects sets of links between a given set of objects and their counterparts.

        This method allows an application to pass in one or more objects of type GCF, Spectrum,
        or MolecularFamily, and retrieve the set of links to their counterpart objects (e.g.
        GCF for Spectrum input) that exceed some threshold score. Returned scores can then
        be accessed using NPLinker.all_links or NPLinker.links_for_obj. 

        Args:
            objects: An instance of one of the above types, or a list containing any mix of them.
            datalinks: The DataLinks instance used for scoring. Defaults to None, in which case the 
                standard internal DataLinks instance is used. 
            scoring_method (str): selects the scoring method to be used (e.g. 'metcalf', 'likescore')
            sig_percentile (int): value between 0-100 indicating the percentile value calculated from the
                randomised scoring matrix that candidate links must exceed to be returned. May be None,
                in which case you must provide a value for scoring_cutoff.
            scoring_cutoff (float): an absolute threshold that candidate links must exceed to be
                returned. Must be None if sig_percentile is not None (i.e. the two are mutually 
                exclusive). 

        Returns:
            If either a single instance or a list with uniform type is used as input, the return value
            is a list containing the objects for which links were found to satisfy the scoring criteria.
            In the case where a mixed list of objects is used, the return value will instead be a 3 element
            list, containing the individual results for the classes GCF, Spectrum, and MolecularFamily 
            (in that order). 
        """
        # single object, pass it straight through
        if not isinstance(objects, list):
            logger.debug('get_links: single object')
            return self._get_links([objects], datalinks, scoring_method, sig_percentile, scoring_cutoff)

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
            return self._get_links(objects, datalinks, scoring_method, sig_percentile, scoring_cutoff)

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
                results[i].append(self._get_links(obj_lists[i], datalinks, scoring_method, sig_percentile, scoring_cutoff))

        return results

    def _store_object_links(self, src, dst, score, scoring_method):
        # structure of self._links is:
        # src objects as keys, values are scoring methods
        #   scoring methods as keys, values are dst objects
        #       dst objects as keys, values are link scores
        if src not in self._links:
            self._links[src] = {m: {} for m in NPLinker.SCORING_METHODS}

        if scoring_method not in self._links[src]: 
            raise Exception('Should never happen!')

        self._links[src][scoring_method][dst] = score

    def _get_links(self, objects, datalinks, scoring_method, sig_percentile, scoring_cutoff):
        if len(objects) == 0:
            return []

        if sig_percentile is not None and scoring_cutoff is not None:
            raise Exception('Cannot use both sig_percentile and scoring_cutoff, one of the two must be None')

        if self._linkfinder is None:
            raise Exception('Need to call process_dataset first')

        if scoring_method not in NPLinker.SCORING_METHODS:
            raise Exception('unknown scoring method "{}"'.format(scoring_method))

        datalinks = datalinks or self._datalinks
        input_type = type(objects[0])

        # if using the randomised scores to generate a significance threshold for the links, 
        # then just call get_links without a cutoff value so it will return all scores
        # for the given objects, which we can then apply the perecentile threshold to (this
        # seemed easier than modifying the LinkFinder implementation)
        if sig_percentile is not None:
            scoring_cutoff = None
            if sig_percentile < 0 or sig_percentile > 100:
                raise Exception('sig_percentile invalid! Expected 0-100, got {}'.format(sig_percentile))

        results = self._linkfinder.get_links(datalinks, objects, scoring_method, scoring_cutoff)

        if sig_percentile is not None:
            # if using percentile option, need to filter the links based on that before continuing...
            results = self._get_links_percentile(input_type, objects, results, scoring_method, sig_percentile)

        scores_found = set()

        if input_type == GCF:
            logger.debug('get_links: result type is GCF, inputs={}, results={}'.format(len(objects), len(results)))
            # for GCF input, results contains two arrays of shape (3, x), 
            # which contain spec-gcf and fam-gcf links respectively 
            result_gcf_spec, result_gcf_fam = results[0], results[1]

            for res, type_ in [(result_gcf_spec, Spectrum), (result_gcf_fam, MolecularFamily)]:
                if res.shape[1] == 0:
                    logger.debug('Found no links for {}, [{}]'.format(objects, type_))
                    continue # no results


                # TODO: self._families is never populated atm, so this breaks with MolecularFamily results!
                if type_ == MolecularFamily:
                    logger.warning('FIXME: MolecularFamily support')
                    continue
                logger.debug('Found {} links for {}, [{}]'.format(res.shape[1], objects, type_))

                lookup = {gcf.id: objects.index(gcf) for gcf in objects}

                # for each entry in the results (each Spectrum or MolecularFamily)
                for j in range(res.shape[1]):
                    # extract the ID of the object and get the object itself
                    obj_id = int(res[NPLinker.R_DST_ID, j])
                    obj = self._spectra[obj_id] if type_ == Spectrum else self._families[obj_id]

                    # retrieve the GCF object too
                    index = lookup[res[NPLinker.R_SRC_ID, j]]
                    gcf = objects[index]

                    # finally, create the entry in self._links and record that this GCF
                    # has at least one link associated with it
                    self._store_object_links(gcf, obj, res[NPLinker.R_SCORE, j], scoring_method)
                    scores_found.add(gcf)

        else:
            logger.debug('get_links: result type is Spectrum/MolFam, inputs={}, results={}'.format(len(objects), len(results)))
            # for non-GCF input, result is a list containing a single array, shape (3, x)
            # where x is the total number of links found
            results = results[0]
            if results.shape[1] == 0:
                logger.debug('Found no links for {}'.format(objects))
                return []

            logger.debug('Found {} links for {}'.format(results.shape[1], objects))

            lookup = {obj.id: objects.index(obj) for obj in objects}

            # for each entry in the results (each GCF)
            for j in range(results.shape[1]):
                # extract the ID of the GCF and use that to get the object
                gcf_id = int(results[NPLinker.R_DST_ID, j])
                gcf = self._gcfs[gcf_id]

                # retrieve the source object similarly
                index = lookup[results[NPLinker.R_SRC_ID, j]]
                obj = objects[index]

                # finally, create the entry in self._links and record that this Spectrum or
                # MolecularFamily has at least one link associated with it
                self._store_object_links(obj, gcf, results[NPLinker.R_SCORE, j], scoring_method)
                scores_found.add(obj)

        return list(scores_found)

    def _get_links_percentile(self, input_type, objects, results, scoring_method, sig_percentile):
        # self._random_scores contains the randomised and stacked scores for each given type of 
        # object, so pull those out and apply np.percentile with the given threshold to the correct 
        # matrices here. the results object from get_links will contain *ALL* potential candiate
        # links at this stage - e.g. for a Spectrum input, it will return all GCFs

        # retrieve the IDS for all the supplied objects
        # TODO will need to add .id to MolecularFamily or this will break
        obj_ids = [obj.id for obj in objects]

        if input_type == GCF:
            # 2-element list returned with spec-gcf, fam-gcf results
            results_spec, results_fam = results

            # lookup the appropriate random scores
            rand_spec_scores = self._random_scores[scoring_method][input_type][Spectrum]
            rand_fam_scores = self._random_scores[scoring_method][input_type][MolecularFamily]

            # calculate the threshold values for each input object
            spec_perc_thresholds = [np.percentile(rand_spec_scores[:, id], sig_percentile) for id in obj_ids]
            fam_perc_thresholds = [np.percentile(rand_fam_scores[:, id], sig_percentile) for id in obj_ids]

            perc_results = [np.zeros((3, 0)), np.zeros((3, 0))]

            data = [(results_spec, spec_perc_thresholds, len(self._spectra))]
            data.append((results_fam, fam_perc_thresholds, len(self._families)))

            for d, (full_results, thresholds, num_objs) in enumerate(data):
                for i in range(len(objects)):
                    # TODO is this the best way to do this??
                    # want to get the indices of full_results where a) the source ID matches the
                    # current object and b) the score exceeds the threshold for that object
                    obj_indices = np.intersect1d(np.where(full_results[NPLinker.R_SCORE, :] >= thresholds[i]),
                                                 np.where(full_results[NPLinker.R_SRC_ID, :] == objects[i].id))

                    # concat the selected columns to the final result
                    perc_results[d] = np.c_[perc_results[d], full_results[:, obj_indices]]

            return [perc_results[0], perc_results[1]]

        elif input_type == Spectrum or input_type == MolecularFamily:
            # here results will be a list which in turn contains a single array, so just extract that
            results = results[0]

            # lookup the appropriate random scores
            rand_scores = self._random_scores[scoring_method][input_type]

            # calculate the threshold values for each input object
            perc_thresholds = [np.percentile(rand_scores[id, :], sig_percentile) for id in obj_ids]

            perc_results = np.zeros((3, 0))

            for i in range(len(objects)):
                # TODO is this the best way to do this??
                # want to get the indices of full_results where a) the source ID matches the
                # current object and b) the score exceeds the threshold for that object
                obj_indices = np.intersect1d(np.where(results[NPLinker.R_SCORE, :] >= perc_thresholds[i]),
                                             np.where(results[NPLinker.R_SRC_ID, :] == obj_ids[i]))

                # build up a new matrix containing only the sufficiently high scoring columns
                perc_results = np.c_[perc_results, results[:, obj_indices]]

            return [perc_results]

        logger.error('Bad input type in _get_links_percentile??')
        return None

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

    @property
    def all_links(self):
        """
        Returns the links for all objects and scoring methods
        """
        return self._links

    def links_for_obj(self, obj, type_=None, scoring_method=DEFAULT_SCORING):
        if type_ is None:
            return [(obj, score) for obj, score in self._links[obj][scoring_method].items()]

        # this is only really useful for GCFs, to pick between MolFam/Spectrum results
        return [(obj, score) for obj, score in self._links[obj][scoring_method].items() if isinstance(obj, type_)]

if __name__ == "__main__":
    # TODO will need to handle arguments here (via argparse)
    # TODO possible to auto-discover some/all of these filenames for a given dataset?

    if len(sys.argv) < 2:
        print('Usage: nplinker.py <dataset root>')
        sys.exit(0)

    # crusemann dataset
    DATASET = sys.argv[1]
    logger.info('Loading dataset from path: {}'.format(DATASET))
        
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

    # generate all datalinks objects etc
    if not npl.process_dataset(random_count=3):
        print('Failed to process dataset')
        sys.exit(-1)

    # get_links with a single Spectrum/MolecularFamily
    obj = npl.spectra[2763:2765]
    result = npl.get_links(obj, scoring_method='metcalf', sig_percentile=99, scoring_cutoff=None)
    # result = npl.get_links(obj, scoring_method='metcalf', sig_percentile=None, scoring_cutoff=135)
    # if the input Spectrum had any links above the cutoff, it will appear in result:
    if len(result) > 0:
        print('Found links to {} objects: {}'.format(len(result), result))

        # print info about GCFs each Spectrum is linked to
        for spec in result:
            gcfs = npl.all_links[spec]['metcalf']
            # or
            gcfs = npl.links_for_obj(spec)
            print('{} has links to {} GCFs'.format(spec, len(gcfs)))

    # get_links with a list of GCFs
    objects = [npl.gcfs[763]] 
    objects = npl.gcfs[763:765]
    result = npl.get_links(objects, scoring_method='metcalf', sig_percentile=99)
    # result = npl.get_links(objects, scoring_method='metcalf', sig_percentile=None, scoring_cutoff=150)
    if len(result) > 0:
        print('{}/{} objects with links found!'.format(len(result), len(objects)))
        for gcf in result:
            spectra = npl.links_for_obj(gcf, Spectrum)
            print('GCF {} has links to {} spectra'.format(gcf, len(spectra)))

    sys.exit(0)

    # finding common_strains
    common_strains = npl.get_common_strains(npl.gcfs[763:766], npl.spectra[3403])
    for spec, gcfs in common_strains.items():
        for gcf, strains in gcfs.items():
            for strain in strains:
                print('{} and {} share strain #{}'.format(spec, gcf, strain))

