import os
import numpy as np
from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.logconfig import LogConfig
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.pickler import load_pickled_data
from nplinker.pickler import save_pickled_data
from nplinker.scoring.methods import ScoringMethod
from nplinker.scoring.object_link import ObjectLink
from nplinker.scoring.linking.data_linking import DataLinks
from nplinker.scoring.linking.link_finder import LinkFinder

logger = LogConfig.getLogger(__file__)


class MetcalfScoring(ScoringMethod):

    DATALINKS = None
    LINKFINDER = None
    NAME = 'metcalf'

    # enumeration for accessing results of LinkFinder.get_links, which are (3, num_links) arrays:
    # - R_SRC_ID: the ID of an object that was supplied as input to get_links
    # - R_DST_ID: the ID of an object that was discovered to have a link to an input object
    # - R_SCORE: the score for the link between a pair of objects
    R_SRC_ID, R_DST_ID, R_SCORE = range(3)

    def __init__(self, npl):
        super().__init__(npl)
        self.cutoff = 1.0
        self.standardised = True

    @staticmethod
    def setup(npl):
        logger.info(
            'MetcalfScoring.setup (bgcs={}, gcfs={}, spectra={}, molfams={}, strains={})'
            .format(len(npl.bgcs), len(npl.gcfs), len(npl.spectra),
                    len(npl.molfams), len(npl.strains)))

        cache_dir = os.path.join(npl.root_dir, 'metcalf')
        cache_file = os.path.join(cache_dir, 'metcalf_scores.pckl')
        os.makedirs(cache_dir, exist_ok=True)

        # the metcalf preprocessing can take a long time for large datasets, so it's
        # better to cache as the data won't change unless the number of objects does

        dataset_counts = [
            len(npl.bgcs),
            len(npl.gcfs),
            len(npl.spectra),
            len(npl.molfams),
            len(npl.strains)
        ]
        datalinks, linkfinder = None, None
        if os.path.exists(cache_file):
            logger.debug('MetcalfScoring.setup loading cached data')
            cache_data = load_pickled_data(npl, cache_file)
            cache_ok = True
            if cache_data is not None:
                (counts, datalinks, linkfinder) = cache_data
                # need to invalidate this if dataset appears to have changed
                for i in range(len(counts)):
                    if counts[i] != dataset_counts[i]:
                        logger.info(
                            'MetcalfScoring.setup invalidating cached data!')
                        cache_ok = False
                        break

            if cache_ok:
                MetcalfScoring.DATALINKS = datalinks
                MetcalfScoring.LINKFINDER = linkfinder

        if MetcalfScoring.DATALINKS is None:
            logger.info(
                'MetcalfScoring.setup preprocessing dataset (this may take some time)'
            )
            MetcalfScoring.DATALINKS = DataLinks()
            MetcalfScoring.DATALINKS.load_data(npl._spectra, npl._gcfs,
                                               npl._strains)
            # TODO fix crash with this set to True, see https://github.com/sdrogers/nplinker/issues/57
            MetcalfScoring.DATALINKS.find_correlations(
                include_singletons=False)
            MetcalfScoring.LINKFINDER = LinkFinder()
            MetcalfScoring.LINKFINDER.metcalf_scoring(MetcalfScoring.DATALINKS,
                                                      type='spec-gcf')
            MetcalfScoring.LINKFINDER.metcalf_scoring(MetcalfScoring.DATALINKS,
                                                      type='fam-gcf')
            logger.debug('MetcalfScoring.setup caching results')
            save_pickled_data((dataset_counts, MetcalfScoring.DATALINKS,
                               MetcalfScoring.LINKFINDER), cache_file)

        logger.info('MetcalfScoring.setup completed')

    @property
    def datalinks(self):
        return MetcalfScoring.DATALINKS

    def _metcalf_postprocess_met(self, linkfinder, results, input_type):
        logger.debug(
            'Postprocessing results for standardised Metcalf scores (met input)'
        )
        # results will be links from EITHER Spectrum OR MolFam => GCF here

        # need to know if the metabolomic objects given as input are Spectrum/MolFam
        met_objs = self.npl.spectra if input_type == Spectrum else self.npl.molfams
        new_src, new_dst, new_sco = [], [], []

        # go through each pair of input objects and calculate their standardised scores
        for i in range(len(results[0][self.R_SRC_ID])):
            met_obj = met_objs[int(results[0][self.R_SRC_ID][i])]
            # met_obj will now be either a Spectrum or a MolecularFamily, but
            # doesn't matter which (in this implementation at least) because they
            # both have a .strains attribute which is the only thing we need. For
            # Spectra it's the number of strains, for a MolFam it's the total
            # number of *unique* strains across all Spectra in that family.
            met_strains = len(met_obj.strains)
            gcf = self.npl.gcfs[int(results[0][self.R_DST_ID][i])]
            gen_strains = len(gcf.strains)

            # lookup expected + variance values based on strain counts
            expected = linkfinder.metcalf_expected[met_strains][gen_strains]
            variance_sqrt = linkfinder.metcalf_variance_sqrt[met_strains][
                gen_strains]

            # calculate the final score based on the basic Metcalf score for these two
            # particular objects
            final_score = (results[0][self.R_SCORE][i] -
                           expected) / variance_sqrt

            # finally apply the scoring cutoff and store the result
            if self.cutoff is None or (final_score >= self.cutoff):
                new_src.append(int(results[0][self.R_SRC_ID][i]))
                new_dst.append(int(results[0][self.R_DST_ID][i]))
                new_sco.append(final_score)

        # overwrite original "results" with equivalent new data structure
        return [np.array([new_src, new_dst, new_sco])]

    def _metcalf_postprocess_gen(self, linkfinder, results, input_type):
        logger.debug(
            'Postprocessing results for standardised Metcalf scores (gen input)'
        )
        # results will be links from GCF to BOTH Spectrum and MolFams here (first
        # element Spectra, second MolFams)

        new_results = []
        met_objs_list = [self.npl.spectra, self.npl.molfams]

        # iterate over the Spectrum results and then the MolFam results
        for m, met_objs in enumerate(met_objs_list):
            new_src, new_dst, new_sco = [], [], []

            # go through each pair of input objects and calculate their standardised scores
            for i in range(len(results[m][self.R_SRC_ID])):
                gcf = self.npl.gcfs[int(results[m][self.R_SRC_ID][i])]
                gen_strains = len(gcf.strains)

                # met_obj will now be either a Spectrum or a MolecularFamily, but
                # doesn't matter which (in this implementation at least) because they
                # both have a .strains attribute which is the only thing we need. For
                # Spectra it's the number of strains, for a MolFam it's the total
                # number of *unique* strains across all Spectra in that family.
                met_obj = met_objs[int(results[m][self.R_DST_ID][i])]
                met_strains = len(met_obj.strains)

                # lookup expected + variance values based on strain counts
                expected = linkfinder.metcalf_expected[met_strains][
                    gen_strains]
                variance_sqrt = linkfinder.metcalf_variance_sqrt[met_strains][
                    gen_strains]

                # calculate the final score based on the basic Metcalf score for these two
                # particular objects
                final_score = (results[m][self.R_SCORE][i] -
                               expected) / variance_sqrt

                # finally apply the scoring cutoff and store the result
                if self.cutoff is None or (final_score >= self.cutoff):
                    new_src.append(int(results[m][self.R_SRC_ID][i]))
                    new_dst.append(int(results[m][self.R_DST_ID][i]))
                    new_sco.append(final_score)

            # overwrite original "results" with equivalent new data structure
            new_results.append(np.array([new_src, new_dst, new_sco]))

        return new_results

    def get_links(self, objects, link_collection):
        # enforce constraint that the list must contain a set of identically typed objects
        if not all(isinstance(x, type(objects[0])) for x in objects):
            raise Exception(
                'MetcalfScoring: uniformly-typed list of objects is required')

        # also can't handle BGCs here, must be one of the other 3 types (GCF/Spectrum/MolecularFamily)
        if isinstance(objects[0], BGC):
            raise Exception(
                'MetcalfScoring requires input type GCF/Spectrum/MolecularFamily, not BGC'
            )

        datalinks = MetcalfScoring.DATALINKS
        linkfinder = MetcalfScoring.LINKFINDER
        input_type = type(objects[0])

        logger.debug('MetcalfScoring: standardised = {}'.format(
            self.standardised))
        if not self.standardised:
            results = linkfinder.get_links(datalinks, objects, self.name,
                                           self.cutoff)
        else:
            # get the basic Metcalf scores BUT ignore the cutoff value here by setting
            # it to None. The actual user-supplied cutoff value is applied further down
            # once the standardised scores for these results have been calculated.
            results = linkfinder.get_links(datalinks, objects, self.name, None)

            # The "results" object varies slightly depending on the input provided
            # to the LinkFinder class:
            #  - given Spectra/MolFam input, it will be a single element list containing
            #   a (3, x) array, where the first row contains source (input) object
            #   IDs, the second contains destination (linked) object IDs, and the
            #   third contains regular Metcalf scores for those pairs of objects.
            #  - however for GCF input, "results" is instead a 2-element list where
            #   each entry has the same structure as described above, with the first
            #   entry describing GCF-Spectrum links and the second GCF-MolFam links.

            gcf_input = (input_type == GCF)

            if not gcf_input:
                results = self._metcalf_postprocess_met(
                    linkfinder, results, input_type)
            else:
                results = self._metcalf_postprocess_gen(
                    linkfinder, results, input_type)

        scores_found = set()
        metcalf_results = {}

        if input_type == GCF:
            logger.debug(
                'MetcalfScoring: input_type=GCF, result_type=Spec/MolFam, inputs={}, results={}'
                .format(len(objects), results[0].shape))
            # for GCF input, results contains two arrays of shape (3, x),
            # which contain spec-gcf and fam-gcf links respectively
            result_gcf_spec, result_gcf_fam = results[0], results[1]

            for res, type_ in [(result_gcf_spec, Spectrum),
                               (result_gcf_fam, MolecularFamily)]:
                if res.shape[1] == 0:
                    if type_ != MolecularFamily:
                        logger.debug(
                            'Found no links for {} input objects (type {})'.
                            format(len(objects), type_))
                    continue  # no results

                # for each entry in the results (each Spectrum or MolecularFamily)
                for j in range(res.shape[1]):
                    # extract the ID of the object and get the object itself
                    obj_id = int(res[self.R_DST_ID, j])
                    obj = self.npl._spectra[
                        obj_id] if type_ == Spectrum else self.npl._molfams[
                            obj_id]

                    # retrieve the GCF object too (can use its internal ID to index
                    # directly into the .gcfs list)
                    gcf = self.npl._gcfs[int(res[self.R_SRC_ID][j])]

                    # record that this GCF has at least one link associated with it
                    scores_found.add(gcf)

                    # save the scores
                    if gcf not in metcalf_results:
                        metcalf_results[gcf] = {}
                    metcalf_results[gcf][obj] = ObjectLink(
                        gcf, obj, self, res[self.R_SCORE, j])

        else:
            logger.debug(
                'MetcalfScoring: input_type=Spec/MolFam, result_type=GCF, inputs={}, results={}'
                .format(len(objects), results[0].shape))
            # for non-GCF input, result is a list containing a single array, shape (3, x)
            # where x is the total number of links found
            results = results[0]
            if results.shape[1] == 0:
                logger.debug('Found no links for {} input objects'.format(
                    len(objects)))
                link_collection._add_links_from_method(self, metcalf_results)
                # can just bail out here in this case
                logger.debug('MetcalfScoring: completed')
                return link_collection

            # for each entry in the results (each GCF)
            for j in range(results.shape[1]):
                # extract the ID of the GCF and use that to get the object itself
                gcf = self.npl._gcfs[int(results[self.R_DST_ID, j])]

                # retrieve the Spec/MolFam object too (can use its internal ID to index
                # directly into the appropriate list)
                obj_id = int(results[self.R_SRC_ID, j])
                obj = self.npl._spectra[
                    obj_id] if input_type == Spectrum else self.npl._molfams[
                        obj_id]

                # record that this Spectrum or MolecularFamily has at least one link associated with it
                scores_found.add(obj)

                # save the scores
                if obj not in metcalf_results:
                    metcalf_results[obj] = {}
                metcalf_results[obj][gcf] = ObjectLink(
                    obj, gcf, self, results[self.R_SCORE, j])

        logger.debug('MetcalfScoring found {} results'.format(
            len(metcalf_results)))
        link_collection._add_links_from_method(self, metcalf_results)
        logger.debug('MetcalfScoring: completed')
        return link_collection

    def format_data(self, data):
        # for metcalf the data will just be a floating point value (i.e. the score)
        return f'{data:.4f}'

    def sort(self, objects, reverse=True):
        # sort based on score
        return sorted(objects,
                      key=lambda objlink: objlink[self],
                      reverse=reverse)
