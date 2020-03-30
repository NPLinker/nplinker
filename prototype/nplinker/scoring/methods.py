import random
from types import SimpleNamespace

import numpy as np

from .data_linking import DataLinks, LinkFinder
from ..genomics import BGC, GCF
from ..metabolomics import Spectrum, MolecularFamily

from ..logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class ScoringMethod(object):

    NAME = 'ScoringMethod'

    def __init__(self, npl):
        self.npl = npl
        self.name = self.__class__.NAME

    @staticmethod
    def setup(npl):
        """Perform any one-off initialisation required (will only be called once)"""
        pass

    def get_links(self, objects):
        """Given a set of objects, return link information"""
        return {}

class TestScoring(ScoringMethod):

    NAME = 'testscore'

    def __init__(self, npl):
        super(TestScoring, self).__init__(npl)
        self.foo = 123

    @staticmethod
    def setup(npl):
        logger.info('TestScoring setup')

    def get_links(self, objects):
        # randomly pick some objects and scores to return
        results = {}
        for obj in objects:
            if random.random() > 0.5:
                results[obj] = (random.random() * 10, self.foo)

        return results

class MetcalfScoring(ScoringMethod):

    DATALINKS = None
    LINKFINDER = None
    NAME = 'metcalf'

    def __init__(self, npl):
        super(MetcalfScoring, self).__init__(npl)
        self.cutoff = 1.0
        self.standardised = True

    @staticmethod
    def setup(npl):
        logger.info('MetcalfScoring.setup')
        MetcalfScoring.DATALINKS = DataLinks()
        MetcalfScoring.DATALINKS.load_data(npl._spectra, npl._gcfs, npl._strains)
        MetcalfScoring.DATALINKS.find_correlations()
        MetcalfScoring.LINKFINDER = LinkFinder()
        MetcalfScoring.LINKFINDER.metcalf_scoring(MetcalfScoring.DATALINKS, type='spec-gcf')
        MetcalfScoring.LINKFINDER.metcalf_scoring(MetcalfScoring.DATALINKS, type='fam-gcf')
        logger.info('MetcalfScoring.setup completed')

    @property
    def datalinks(self):
        return MetcalfScoring.DATALINKS

    def get_links(self, objects):
        # enforce constraint that the list must contain a set of identically typed objects
        if not all(isinstance(x, type(objects[0])) for x in objects):
            raise Exception('MetcalfScoring: uniformly-typed list of objects is required')

        # also can't handle BGCs here, must be one of the other 3 types (GCF/Spectrum/MolecularFamily)
        if isinstance(objects[0], BGC):
            raise Exception('MetcalfScoring requires input type GCF/Spectrum/MolecularFamily, not BGC')

        datalinks = MetcalfScoring.DATALINKS
        linkfinder = MetcalfScoring.LINKFINDER
        input_type = type(objects[0])

        logger.debug('MetcalfScoring: standardised = {}'.format(self.standardised))
        if not self.standardised:
            results = linkfinder.get_links(datalinks, objects, self.name, self.cutoff)
        else:
             # get the basic Metcalf scores BUT ignore the cutoff value here as it should only be applied to 
            # the final scores, not the original ones
            results = linkfinder.get_links(datalinks, objects, self.name, None)

            # TODO molfam support still missing here (as in data_linking.py)

            # results is a (3, x) array for spec/molfam input, where (1, :) gives src obj ID, (2, :) gives
            # dst obj ID, and (3, :) gives scores
            # for GCFs you get [spec, molfam] with above substructure

            # NOTE: need to use sets here because if the "else" cases are executed 
            # can get lots of duplicate object IDs which messes everything up

            spectra_input = (input_type == Spectrum)

            # make type1_objects always == spectra. if input is spectra, just use "objects",
            # otherwise build a list using the object IDs in the results array
            type1_objects = objects if spectra_input else {self.npl.spectra[int(index)] for index in results[0][self.npl.R_DST_ID]}

            # other objs should be "objects" if input is GCF, otherwise create 
            # from the results array as above
            other_objs = objects if not spectra_input else {self.npl._gcfs[int(index)] for index in results[0][self.npl.R_DST_ID]}

            # build a lookup table for pairs of object IDs and their index in the results array
            if spectra_input:
                pairs_lookup = {(results[0][self.npl.R_SRC_ID][i], results[0][self.npl.R_DST_ID][i]): i for i in range(len(results[0][self.npl.R_SRC_ID]))}
            else:
                pairs_lookup = {(results[0][self.npl.R_DST_ID][i], results[0][self.npl.R_SRC_ID][i]): i for i in range(len(results[0][self.npl.R_SRC_ID]))}

            # apply score update to each spec:gcf pairing and update the corresponding entry in results
            # (the metcalf_expected and metcalf_variance matrices should have been calculated already)
            new_src, new_dst, new_res = [], [], []

            # iterating over spectra
            for i, type1_obj in enumerate(type1_objects):
                met_strains = len(type1_obj.dataset_strains)
                # iterating over GCFs
                for j, other_obj in enumerate(other_objs):
                    gen_strains = len(other_obj.dataset_strains) 

                    # lookup expected + variance values based on strain counts 
                    expected = linkfinder.metcalf_expected[met_strains][gen_strains]
                    variance_sqrt = linkfinder.metcalf_variance_sqrt[met_strains][gen_strains]

                    # now need to extract the original score for this particular pair of objects. 
                    # this requires figuring out the correct index into the results array for the pair,
                    # using the lookup table constructed above
                    k = pairs_lookup[(type1_obj.id, other_obj.id)]

                    # calculate the final score based on the basic Metcalf score for these two
                    # particular objects
                    final_score = (results[0][self.npl.R_SCORE][k] - expected) / variance_sqrt

                    # finally apply the scoring cutoff and store the result
                    if self.cutoff is None or (final_score >= self.cutoff):
                        new_src.append(results[0][self.npl.R_SRC_ID][k])
                        new_dst.append(results[0][self.npl.R_DST_ID][k])
                        new_res.append(final_score)
                
            # overwrite original "results" with equivalent new data structure
            results = [np.array([new_src, new_dst, new_res])]

            if input_type == GCF:
                # TODO molfam...
                results.append(np.zeros((3, 0)))

        scores_found = set()
        metcalf_results = {}

        if input_type == GCF:
            logger.debug('MetcalfScoring: input_type=GCF, result_type=Spec/MolFam, inputs={}, results={}'.format(len(objects), results[0].shape))
            # for GCF input, results contains two arrays of shape (3, x), 
            # which contain spec-gcf and fam-gcf links respectively 
            result_gcf_spec, result_gcf_fam = results[0], results[1]

            for res, type_ in [(result_gcf_spec, Spectrum), (result_gcf_fam, MolecularFamily)]:
                if res.shape[1] == 0:
                    logger.debug('Found no links for {} input objects (type {})'.format(len(objects), type_))
                    continue # no results

                # for each entry in the results (each Spectrum or MolecularFamily)
                for j in range(res.shape[1]):
                    # extract the ID of the object and get the object itself
                    obj_id = int(res[self.npl.R_DST_ID, j])
                    obj = self.npl._spectra[obj_id] if type_ == Spectrum else self.npl._molfams[obj_id]

                    # retrieve the GCF object too (can use its internal ID to index
                    # directly into the .gcfs list)
                    gcf = self.npl._gcfs[int(res[self.npl.R_SRC_ID][j])]

                    # record that this GCF has at least one link associated with it
                    scores_found.add(gcf)

                    # save the scores
                    if gcf not in metcalf_results:
                        metcalf_results[gcf] = []
                    metcalf_results[gcf].append(SimpleNamespace(**{'src': gcf, 'dst': obj, 'score': res[self.npl.R_SCORE, j]}))

        else:
            logger.debug('MetcalfScoring: input_type=Spec/MolFam, result_type=GCF, inputs={}, results={}'.format(len(objects), results[0].shape))
            # for non-GCF input, result is a list containing a single array, shape (3, x)
            # where x is the total number of links found
            results = results[0]
            if results.shape[1] == 0:
                logger.debug('Found no links for {} input objects'.format(len(objects)))
                return []

            # for each entry in the results (each GCF)
            for j in range(results.shape[1]):
                # extract the ID of the GCF and use that to get the object itself
                gcf = self.npl._gcfs[int(results[self.npl.R_DST_ID, j])]

                # retrieve the Spec/MolFam object too (can use its internal ID to index
                # directly into the appropriate list)
                obj_id = int(results[self.npl.R_SRC_ID, j])
                obj = self.npl._spectra[obj_id] if input_type == Spectrum else self.npl._molfams[obj_id]

                # record that this Spectrum or MolecularFamily has at least one link associated with it
                scores_found.add(obj)

                # save the scores
                if obj not in metcalf_results:
                    metcalf_results[obj] = []
                metcalf_results[obj].append(SimpleNamespace(**{'src': obj, 'dst': gcf, 'score': results[self.npl.R_SCORE, j]}))

        return metcalf_results
