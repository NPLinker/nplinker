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

import itertools
import os
import random
import numpy as np
from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.logconfig import LogConfig
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from nplinker.pickler import load_pickled_data
from nplinker.pickler import save_pickled_data
from nplinker.scoring.object_link import ObjectLink
from nplinker.scoring.rosetta.rosetta import Rosetta
from .linking.data_linking import DataLinks
from .linking.data_linking import LinkFinder


logger = LogConfig.getLogger(__file__)

# TODO update/expand comments in this file!

class ScoringMethod():

    NAME = 'ScoringMethod'

    def __init__(self, npl):
        self.npl = npl
        self.name = self.__class__.NAME

    @staticmethod
    def setup(npl):
        """Perform any one-off initialisation required (will only be called once)"""
        pass

    def get_links(self, objects, link_collection):
        """Given a set of objects, return link information"""
        return link_collection

    def format_data(self, data):
        """Given whatever output data the method produces, return a readable string version"""
        return ''

    def sort(self, objects, reverse=True):
        """Given a list of objects, return them sorted by link score"""
        return objects


class TestScoring(ScoringMethod):

    NAME = 'testscore'

    def __init__(self, npl):
        super().__init__(npl)
        self.value = 0.5
        self.mc = MetcalfScoring(npl)

    @staticmethod
    def setup(npl):
        logger.info('TestScoring setup')

    def get_links(self, objects, link_collection):
        mc_results = self.mc.get_links(objects, link_collection)
        num_to_keep = int(len(mc_results) * self.value)
        results = {
            obj: data
            for obj, data in list(mc_results.links.items())[:num_to_keep]
        }
        for links in results.values():
            for link in links.values():
                # this is just to make things work properly for the test
                # method, shouldn't do stuff like this normally
                link._method_data[self] = random.random()
                del link._method_data[self.mc]

        logger.debug(f'TestScoring found {len(results)} results')
        link_collection._add_links_from_method(self, results)
        return link_collection

    def format_data(self, data):
        return self.mc.format_data(data)

    def sort(self, objects, reverse=True):
        # nothing
        return objects


class RosettaScoring(ScoringMethod):

    NAME = 'rosetta'
    ROSETTA_OBJ = None

    def __init__(self, npl):
        super().__init__(npl)
        self.bgc_to_gcf = True

        self.spec_score_cutoff = 0.0
        self.bgc_score_cutoff = 0.0

    @staticmethod
    def setup(npl):
        logger.info('RosettaScoring setup')
        RosettaScoring.ROSETTA_OBJ = Rosetta(npl, ignore_genomic_cache=False)
        ms1_tol = Rosetta.DEF_MS1_TOL
        ms2_tol = Rosetta.DEF_MS2_TOL
        score_thresh = Rosetta.DEF_SCORE_THRESH
        min_match_peaks = Rosetta.DEF_MIN_MATCH_PEAKS

        # allow overridding params via config file
        config = npl.config
        if 'scoring' in config and 'rosetta' in config['scoring']:
            rc = config['scoring']['rosetta']
            ms1_tol = rc.get('ms1_tol', Rosetta.DEF_MS1_TOL)
            ms2_tol = rc.get('ms2_tol', Rosetta.DEF_MS2_TOL)
            score_thresh = rc.get('score_thresh', Rosetta.DEF_SCORE_THRESH)
            min_match_peaks = rc.get('min_match_peaks',
                                     Rosetta.DEF_MIN_MATCH_PEAKS)

        RosettaScoring.ROSETTA_OBJ.run(npl.spectra, npl.bgcs, ms1_tol, ms2_tol,
                                       score_thresh, min_match_peaks)
        logger.info('RosettaScoring setup completed')

    def _include_hit(self, hit):
        if hit.spec_match_score < self.spec_score_cutoff or hit.bgc_match_score < self.bgc_score_cutoff:
            return False

        return True

    def _insert_result_gen(self, results, src, hit):
        if src not in results:
            results[src] = {}
        # Rosetta can produce multiple "hits" per link, need to
        # ensure the ObjectLink contains all the RosettaHit objects
        # in these cases
        if hit.spec in results[src]:
            original_data = results[src][hit.spec].data(self)
            results[src][hit.spec].set_data(self, original_data + [hit])
        else:
            results[src][hit.spec] = ObjectLink(src,
                                                hit.spec,
                                                self,
                                                data=[hit])

        return results

    def _insert_result_met(self, results, spec, target, hit):
        if spec not in results:
            results[spec] = {}
        # Rosetta can produce multiple "hits" per link, need to
        # ensure the ObjectLink contains all the RosettaHit objects
        # in these cases
        if target in results[spec]:
            original_data = results[spec][target].data(self)
            results[spec][target].set_data(self, original_data + [hit])
        else:
            results[spec][target] = ObjectLink(spec, target, self, data=[hit])

        return results

    def get_links(self, objects, link_collection):
        # enforce constraint that the list must contain a set of identically typed objects
        if not all(isinstance(x, type(objects[0])) for x in objects):
            raise Exception(
                'RosettaScoring: uniformly-typed list of objects is required')

        if isinstance(objects[0], MolecularFamily):
            raise Exception(
                'RosettaScoring requires input type Spectrum (found MolecularFamily)'
            )

        if isinstance(objects[0], GCF):
            # assume user wants to use all BGCs from these GCFs
            bgcs = list(
                set(itertools.chain.from_iterable(x.bgcs for x in objects)))
            logger.info(
                'RosettaScoring got {} GCFs input, converted to {} BGCs'.
                format(len(objects), len(bgcs)))
            objects = bgcs

        # list of RosettaHit objects which satisfy the current cutoffs
        ro_hits = list(
            filter(lambda hit: self._include_hit(hit),
                   RosettaScoring.ROSETTA_OBJ._rosetta_hits))

        # TODO this might need to be faster
        results = {}
        if isinstance(objects[0], BGC):
            for bgc in objects:
                for hit in ro_hits:
                    if bgc.id == hit.bgc.id:
                        if not self.bgc_to_gcf:
                            # can use the BGC directly
                            results = self._insert_result_gen(
                                results, bgc, hit)
                        else:
                            # if we want the results to contain GCFs instead, need
                            # to iterate over the set of bgc.parents (most will only
                            # have one, hybrid BGCs will have more)
                            for gcf in bgc.parents:
                                results = self._insert_result_gen(
                                    results, gcf, hit)
        else:  # Spectrum
            for spec in objects:
                for hit in ro_hits:
                    if spec.id == hit.spec.id:
                        if not self.bgc_to_gcf:
                            # can use the BGC directly
                            results = self._insert_result_met(
                                results, spec, hit.bgc, hit)
                        else:
                            # if we want the results to contain GCFs instead, need
                            # to iterate over the set of bgc.parents (most will only
                            # have one, hybrid BGCs will have more)
                            for gcf in hit.bgc.parents:
                                results = self._insert_result_met(
                                    results, spec, gcf, hit)

        link_collection._add_links_from_method(self, results)
        logger.debug(f'RosettaScoring found {len(results)} results')
        return link_collection

    def format_data(self, data):
        # TODO
        return f'{len(data)} hits'

    def sort(self, objects, reverse=True):
        # TODO
        return objects


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


