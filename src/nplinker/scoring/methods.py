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
import time
import numpy as np
from ..genomics import BGC
from ..genomics import GCF
from ..logconfig import LogConfig
from ..metabolomics import MolecularFamily
from ..metabolomics import Spectrum
from ..pickler import load_pickled_data
from ..pickler import save_pickled_data
from ..scoring.rosetta.rosetta import Rosetta
from .data_linking import DataLinks
from .data_linking import LinkFinder


logger = LogConfig.getLogger(__file__)

# TODO update/expand comments in this file!


class LinkCollection():
    """
    Class which stores the results of running one or more scoring methods.

    It provides access to the set of objects which were found to have links,
    the set of objects linked to each of those objects, and the information
    produced by the scoring method(s) about each link.

    There are also some useful utility methods to filter the original results.
    """

    def __init__(self, and_mode=True):
        self._methods = set()
        self._link_data = {}
        self._targets = {}
        self._and_mode = and_mode

    def _add_links_from_method(self, method, object_links):
        if method in self._methods:
            # this is probably an error...
            raise Exception(
                'Duplicate method found in LinkCollection: {}'.format(
                    method.name))

        # if this is the first set of results to be generated, can just dump
        # them all straight in
        if len(self._methods) == 0:
            self._link_data = {k: v for k, v in object_links.items()}
        else:
            # if already some results added, in OR mode can just merge the new set
            # with the existing set, but in AND mode need to ensure we end up with
            # only results that appear in both sets

            if not self._and_mode:
                logger.debug(
                    'Merging {} results from method {} in OR mode'.format(
                        len(object_links), method.name))
                self._merge_or_mode(object_links)
            else:
                logger.debug(
                    'Merging {} results from method {} in AND mode'.format(
                        len(object_links), method.name))
                self._merge_and_mode(object_links)

        self._methods.add(method)

    def _merge_and_mode(self, object_links):
        # set of ObjectLinks common to existing + new results
        intersect1 = self._link_data.keys() & object_links.keys()

        # iterate over the existing set of link info, remove entries for objects
        # that aren't common to both that and the new set of info, and merge in
        # any common links
        to_remove = set()
        for source, existing_links in self._link_data.items():
            if source not in intersect1:
                to_remove.add(source)
                continue

            links_to_merge = object_links[source]
            intersect2 = existing_links.keys() & links_to_merge.keys()

            self._link_data[source] = {
                k: v
                for k, v in existing_links.items() if k in intersect2
            }

            for target, object_link in object_links[source].items():
                if target in self._link_data[source]:
                    self._link_data[source][target]._merge(object_link)

            if len(self._link_data[source]) == 0:
                to_remove.add(source)

        for source in to_remove:
            del self._link_data[source]

    def _merge_or_mode(self, object_links):
        # source = GCF/Spectrum, links = {Spectrum/GCF: ObjectLink} dict
        for source, links in object_links.items():

            # update the existing dict with the new entries that don't appear in it already
            if source not in self._link_data:
                self._link_data[source] = links
            else:
                self._link_data[source].update({
                    k: v
                    for k, v in links.items()
                    if k not in self._link_data[source]
                })

            # now merge the remainder (common to both)
            for target, object_link in links.items():
                self._link_data[source][target]._merge(object_link)

    def filter_no_shared_strains(self):
        len_before = len(self._link_data)
        self.filter_links(lambda x: len(x.shared_strains) > 0)
        logger.debug('filter_no_shared_strains: {} => {}'.format(
            len_before, len(self._link_data)))

    def filter_sources(self, callable_obj):
        len_before = len(self._link_data)
        self._link_data = {
            k: v
            for k, v in self._link_data.items() if callable_obj(k)
        }
        logger.debug('filter_sources: {} => {}'.format(len_before,
                                                       len(self._link_data)))

    def filter_targets(self, callable_obj, sources=None):
        to_remove = []
        sources_list = self._link_data.keys() if sources is None else sources
        for source in sources_list:
            self._link_data[source] = {
                k: v
                for k, v in self._link_data[source].items() if callable_obj(k)
            }
            # if there are now no links for this source, remove it completely
            if len(self._link_data[source]) == 0:
                to_remove.append(source)

        for source in to_remove:
            del self._link_data[source]

    def filter_links(self, callable_obj, sources=None):
        to_remove = []
        sources_list = self._link_data.keys() if sources is None else sources
        for source in sources_list:
            self._link_data[source] = {
                k: v
                for k, v in self._link_data[source].items() if callable_obj(v)
            }
            # if there are now no links for this source, remove it completely
            if len(self._link_data[source]) == 0:
                to_remove.append(source)

        for source in to_remove:
            del self._link_data[source]

    def get_sorted_links(self, method, source, reverse=True, strict=False):
        # This method allows for the sorting of a set of links according to the
        # sorting implemented by a specific method. However because there may be
        # links from multiple methods present in the collection, it isn't as simple
        # as running <method>.sort(links) and returning the result, because that
        # will only work on links which have the expected method data. To get around
        # this, the "strict" parameter is used. If set to True, it simply returns
        # the sorted links *for the specific method only*, which may be a subset
        # of the total collection if multiple methods were used to generate it. If
        # set to False, it will return a list consisting of the sorted links for
        # the given method, with any remaining links appended in arbitrary order.

        # run <method>.sort on the links found by that method
        sorted_links_for_method = method.sort([
            link for link in self._link_data[source].values()
            if method in link.methods
        ], reverse)

        if not strict:
            # append any remaining links
            sorted_links_for_method.extend([
                link for link in self._link_data[source].values()
                if method not in link.methods
            ])

        return sorted_links_for_method

    def get_all_targets(self):
        return list(
            set(
                itertools.chain.from_iterable(
                    self._link_data[x].keys()
                    for x in self._link_data.keys())))

    @property
    def methods(self):
        return self._methods

    @property
    def sources(self):
        # the set of objects supplied as input, which have links
        return list(self._link_data.keys())

    @property
    def links(self):
        return self._link_data

    @property
    def source_count(self):
        return len(self._link_data)

    @property
    def method_count(self):
        return len(self._methods)

    def __len__(self):
        return len(self._link_data)


class ObjectLink():
    """
    Class which stores information about a single link between two objects.

    There will be at most one instance of an ObjectLink for a given pair of
    objects (source, target) after running 1 or more scoring methods. Some
    methods, e.g. Metcalf, will always produce a single output per link.
    However other methods like Rosetta may find multiple "hits" for a given
    pair. In either case the data for a given method is associated with the
    ObjectLink so it can be retrieved afterwards.

    The information stored is basically:
     - the "source" of the link (original object provided as part of the input)
     - the "target" of the link (linked object, as determined by the method(s) used)
     - a (possibly empty) list of Strain objects shared between source and target
     - the output of the scoring method(s) used for this link (e.g. a metcalf score)
    """

    def __init__(self, source, target, method, data=None, shared_strains=[]):
        self.source = source
        self.target = target
        self.shared_strains = shared_strains
        self._method_data = {method: data}

    def _merge(self, other_link):
        self._method_data.update(other_link._method_data)
        return self

    def set_data(self, method, newdata):
        self._method_data[method] = newdata

    @property
    def method_count(self):
        return len(self._method_data)

    @property
    def methods(self):
        return list(self._method_data.keys())

    def data(self, method):
        return self._method_data[method]

    def __getitem__(self, name):
        if name in self._method_data:
            return self._method_data[name]

        return object.__getitem__(self, name)

    def __hash__(self):
        # return the nplinker internal ID as hash value (for set/dict etc)
        return self.source.id

    def __str__(self):
        return 'ObjectLink(source={}, target={}, #methods={})'.format(
            self.source, self.target, len(self._method_data))

    def __repr__(self):
        return str(self)


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
        super(TestScoring, self).__init__(npl)
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

        logger.debug('TestScoring found {} results'.format(len(results)))
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
        super(RosettaScoring, self).__init__(npl)
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
        logger.debug('RosettaScoring found {} results'.format(len(results)))
        return link_collection

    def format_data(self, data):
        # TODO
        return '{} hits'.format(len(data))

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
        super(MetcalfScoring, self).__init__(npl)
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
        return '{:.4f}'.format(data)

    def sort(self, objects, reverse=True):
        # sort based on score
        return sorted(objects,
                      key=lambda objlink: objlink[self],
                      reverse=reverse)


class NPClassScoring(ScoringMethod):
    NAME = 'npclassscore'

    def __init__(self, npl):
        super(NPClassScoring, self).__init__(npl)
        self.cutoff = 0.25
        self.method_options = npl.chem_classes.class_predict_options
        self.method = self.method_options[0]
        self.num_results = 1  # how many scores do you want for each link
        # take care what targets are as this method can link both bgc/gcf to
        # both spectrum/gcf
        self.equal_targets = False  # if obj GCF, target is spec not MF
        self.both_targets = False  # if obj GCF, target both spec and MF
        # filter out spectra without a score due to missing spectrum classes
        self.filter_missing_scores = False
        self._target_no_scores = set()

    def _npclass_score(self,
                       obj,
                       target,
                       method='mix',
                       obj_classes=None,
                       target_classes=None):
        """Return sorted class link scores for scoring obj and target

        The input objects can be any mix of the following NPLinker types:
            - BGC
            - GCF
            - Spectrum
            - MolecularFamily

        Classes will be retrieved from NPLinker object with _get_gen_classes
        and _get_met_classes or the classes can be predetermined and given as
        input with the function, this should correspond to obj and
        target (if obj is BGC then  obj_classes should be classes for the BGC).
        Can be:
            - BGC/GCF classes: {'as_classes': list(str), 'bigscape_class': str}
            - Spectrum/Molfam classes: tuple of (classes, class_names_indices),
            where classes is a list of list of tuples/None, where each
            tuple is a class and a score (str, float), and class_names_indices
            is a list of ints that relate to the name of a class ontology lvl

        Args:
            obj: one of the possible input objects
            target: one of the possible input objects
            method: str, which classification method to use for spectra. default is 'main'
                options:
                    -'main': use the main method - currently canopus
                    -'mix': use main method first, when no main classes present, use the others
                    if present
                    -'canopus': use only canopus class predictions
                    -'molnetenhancer': use only molnetenhancer class predictions
            obj_classes: default - None, or classes for one of the input types
            target_classes: default - None, or classes for one of the input types

        Returns:
            List of tuple of
            (score, obj_class_lvl, target_class_lvl, obj_class, target_class)
            (float, str, str, str, str)
            List will be empty if either BGC or spectrum classes are missing
        """
        # todo: make subroutines
        # assess what is obj and target
        spec_like = obj
        bgc_like = target
        spec_like_classes_tup = obj_classes
        bgc_like_classes_dict = target_classes
        bgc_to_spec = False
        if isinstance(obj, BGC) or isinstance(obj, GCF):
            bgc_like = obj
            spec_like = target
            bgc_like_classes_dict = obj_classes
            spec_like_classes_tup = target_classes
            bgc_to_spec = True

        # assess method - move to get_links?
        assert method in self.method_options, \
            (f"NPClass method should be one of method options: {self.method_options}, if your method is not " +
             "in the options check if the class predictions (canopus, etc.) are loaded correctly")

        # gather correct classes if not provided, dict for bgcs, tup for spec
        if not bgc_like_classes_dict:
            bgc_like_classes_dict = self._get_gen_classes(bgc_like)
        if not spec_like_classes_tup:
            spec_like_classes_tup = self._get_met_classes(spec_like, method)
        # unpack spec_like classes - both are lists
        spec_like_classes, spec_like_classes_names_inds = spec_like_classes_tup

        scores = [
        ]  # this will be returned if one of the class sides is absent
        std_score = 0  # if link not recorded in scores (mibig) return this score
        # loop through classes that are possible to link (names in class_match object)
        for bgc_class_name in self.npl.class_matches.bgc_class_names:
            if bgc_class_name == "mibig_classes":
                # treat specially as bigscape class needs to be translated to mibig class
                bigscape_class = bgc_like_classes_dict["bigscape_class"]
                # convert bigscape class to mibig class
                bgc_like_classes = [
                    self.npl.class_matches.bigscape_mibig_conversion.get(
                        bigscape_class)
                ]
            else:
                bgc_like_classes = bgc_like_classes_dict.get(bgc_class_name)
            if bgc_like_classes and spec_like_classes:  # check for classes from both sides
                for bgc_class in bgc_like_classes:
                    for chem_class_name in self.npl.class_matches.chem_class_names:
                        # does info exist for this spectrum class level, return index for class level
                        spec_class_level_i = spec_like_classes_names_inds.get(
                            chem_class_name)
                        if spec_class_level_i:
                            spec_class_tup = spec_like_classes[
                                spec_class_level_i]
                            if spec_class_tup:  # if there is a class at this lvl
                                # is a tuple of (name, score) so take [0]
                                spec_class = spec_class_tup[0]
                                # determine direction of scoring: BGC -> spectrum
                                if bgc_to_spec:
                                    score = self.npl.class_matches.\
                                        class_matches[bgc_class_name]\
                                        [chem_class_name].get(bgc_class, {})\
                                        .get(spec_class, std_score)
                                    result_tuple = (score, bgc_class_name,
                                                    chem_class_name, bgc_class,
                                                    spec_class)
                                else:  # spectrum -> BGC
                                    score = self.npl.class_matches.\
                                        class_matches[chem_class_name]\
                                        [bgc_class_name].get(spec_class, {})\
                                        .get(bgc_class, std_score)
                                    result_tuple = (score, chem_class_name,
                                                    bgc_class_name, spec_class,
                                                    bgc_class)
                                scores.append(result_tuple)
        return sorted(scores, reverse=True)

    def _get_targets(self, test_id):
        """
        Get the targets based upon instance of test_id, returns list of targets

        Args:
            test_id: one of the NPLinker objects: BGC, GCF, Spectrum, Molfam
        Returns:
            List of one or more of one of the NPLinker objects
        """
        if isinstance(test_id, BGC) or isinstance(test_id, GCF):  # genome side
            if self.both_targets:  # no matter BGC or GCF take both spec and MF
                targets = self.npl.spectra + self.npl.molfams
            elif isinstance(test_id, BGC):  # obj are BGC
                if self.equal_targets:  # take
                    targets = self.npl.spectra
                else:
                    targets = self.npl.molfams
            else:  # obj are GCF
                if self.equal_targets:
                    targets = self.npl.molfams
                else:
                    targets = self.npl.spectra
        else:  # metabolome side
            if self.both_targets:
                targets = self.npl.bgcs + self.npl.gcfs
            elif isinstance(test_id, Spectrum):
                if self.equal_targets:
                    targets = self.npl.bgcs
                else:
                    targets = self.npl.gcfs
            else:  # obj are molfam
                if self.equal_targets:
                    targets = self.npl.gcfs
                else:
                    targets = self.npl.bgcs
        return targets

    def _get_gen_classes(self, bgc_like, gcf_as_cutoff=0.5):
        """Get classes for genomics objects

        Args:
            bgc_like: BGC or GCF object from NPLinker input objects
            gcf_as_cutoff: float - if GCF, get antismash classes with
                class_matches.get_gcf_as_classes if the class occurs in more
                than this fraction of the GCF, default = 0.5
        Returns:
            Genome-based classes as a dict with antismash and bigscpae classes
            {'as_classes': list(str), 'bigscape_class': str}
        """
        # assess if bgc or gcf
        is_bgc = isinstance(bgc_like, BGC)
        if is_bgc:
            # get parent gcf for bgc
            bgc_like_gcf = [
                gcf for gcf in self.npl.gcfs
                if bgc_like.id in [b.id for b in gcf.bgcs]
            ][0]
            # gather AS classes and convert to names in scoring dict
            as_classes = self.npl.class_matches.convert_as_classes(
                bgc_like.product_prediction.split("."))
            bgc_like_classes_dict = {
                "bigscape_class": bgc_like_gcf.bigscape_class,
                # str - always one bigscape class right?
                "as_classes": as_classes
            }  # list(str)
        else:
            as_classes = self.npl.class_matches.convert_as_classes(
                self.npl.class_matches.get_gcf_as_classes(
                    bgc_like, gcf_as_cutoff))
            bgc_like_classes_dict = {
                "bigscape_class": bgc_like.bigscape_class,
                # str - always one bigscape class right?
                "as_classes": as_classes
            }  # list(str)
        return bgc_like_classes_dict

    def _get_met_classes(self, spec_like, method='mix'):
        """Get chemical classes for a Spectrum or MolFam based on method

        Args:
            spec_like: Spectrum or MolFam, one of the NPLinker input types
            method: str, one of the appropriate methods for chemical class
                predictions (mix, canopus...), default='mix'
        Returns:
            tuple of (classes, class_names_indices),
            where classes is a list of list of tuples/None, where each
            tuple is a class and a score (str, float), and class_names_indices
            is a list of ints that relate to the name of a class ontology lvl
        """
        # assess if spectrum or molfam
        is_spectrum = isinstance(spec_like, Spectrum)

        # gather classes for spectra, using right method
        # choose the main method here by including it as 'main' in the method parameter
        use_canopus = ('main' in method or 'canopus' in method
                       or 'mix' in method) and 'canopus' in self.method_options
        use_mne = ('molnetenhancer' in method or 'mix' in method) and \
                  'molnetenhancer' in self.method_options
        spec_like_classes, spec_like_classes_names, \
            spec_like_classes_names_inds = (None, None, None)
        # the order in which the classes are read, determines the priority (now: first canopus, then mne)
        if use_canopus and not spec_like_classes:
            if is_spectrum:
                # list of list of tuples/None - todo: add to spectrum object?
                # take only 'best' (first) classification per ontology level
                all_classes = self.npl.chem_classes.canopus. \
                    spectra_classes.get(str(spec_like.spectrum_id))
                if all_classes:
                    spec_like_classes = [
                        cls_per_lvl for lvl in all_classes
                        for i, cls_per_lvl in enumerate(lvl) if i == 0
                    ]
                spec_like_classes_names_inds = self.npl.chem_classes.canopus. \
                    spectra_classes_names_inds
            else:  # molfam
                fam_id = str(spec_like.family_id)
                if fam_id == '-1':  # account for singleton families
                    fam_id += f'_{spec_like.spectra[0].spectrum_id}'
                all_classes = self.npl.chem_classes.canopus.molfam_classes.get(
                    fam_id)
                if all_classes:
                    spec_like_classes = [
                        cls_per_lvl for lvl in all_classes
                        for i, cls_per_lvl in enumerate(lvl) if i == 0
                    ]
                spec_like_classes_names_inds = self.npl.chem_classes.canopus. \
                    molfam_classes_names_inds
        if use_mne and not spec_like_classes:
            # if mne or when main/canopus does not get classes
            if is_spectrum:
                spec_like_classes = self.npl.chem_classes.molnetenhancer. \
                    spectra_classes(spec_like.spectrum_id)
            else:  # molfam
                fam_id = str(spec_like.family_id)
                if fam_id == '-1':  # account for singleton families
                    fam_id += f'_{spec_like.spectra[0].spectrum_id}'
                spec_like_classes = self.npl.chem_classes.molnetenhancer. \
                    molfam_classes.get(fam_id)
            # classes are same for molfam and spectrum so names are irrespective of is_spectrum
            spec_like_classes_names_inds = self.npl.chem_classes.\
                molnetenhancer.spectra_classes_names_inds
        return spec_like_classes, spec_like_classes_names_inds

    @staticmethod
    def setup(npl):
        """Perform any one-off initialisation required (will only be called once)"""
        logger.info("Set up NPClassScore scoring")
        met_options = npl.chem_classes.class_predict_options
        logger.info(f"Please choose one of the methods from {met_options}")
        if not met_options:
            logger.warn(f"There are no methods available! This is probably "
                        f"because no class predictions (canopus, etc.) were "
                        f"present or they are not loaded correctly")
        else:
            logger.info(f"Currently the method '{met_options[0]}' is selected")
        # todo: give info about parameters

    def get_links(self, objects, link_collection):
        """Given a set of objects, return link information"""
        # todo: pickle results
        logger.info("Running NPClassScore...")
        begin = time.time()
        first_obj = objects[0]
        targets = self._get_targets(first_obj)
        obj_is_gen = isinstance(first_obj, BGC) or isinstance(first_obj, GCF)

        # only get target classes once for each target here
        if obj_is_gen:  # obj is genome so get metabolome classes for target
            targets_classes = [
                self._get_met_classes(target, self.method)
                for target in targets
            ]
        else:
            targets_classes = [
                self._get_gen_classes(target) for target in targets
            ]

        logger.info("Using Metcalf scoring to get shared strains")
        # get mapping of shared strains
        if not self.npl._datalinks:
            self.npl._datalinks = self.npl.scoring_method(
                MetcalfScoring.NAME).datalinks
        # this is a dict with structure:
        #   tup(Spectrum/MolecularFamily, GCF) => array of strain indices
        common_strains = self.npl._datalinks.common_strains(
            objects, targets, True)
        logger.info(f"Calculating NPClassScore for {len(objects)} objects to "
                    f"{len(targets)} targets ({len(common_strains)} pairwise "
                    f"interactions that share at least 1 strain). This might "
                    f"take a while.")

        results = {}
        for obj in objects:
            results[obj] = {}
            # get obj class
            if obj_is_gen:
                obj_classes = self._get_gen_classes(obj)
            else:
                obj_classes = self._get_met_classes(obj, self.method)

            for target, target_classes in zip(targets, targets_classes):
                # only consider targets that have shared strains
                common_tup = (obj, target)
                if obj_is_gen:
                    common_tup = (target, obj)
                shared_strains = common_strains.get(common_tup)
                if shared_strains is not None:  # todo: fix for bgcs
                    full_score = self._npclass_score(
                        obj, target, self.method, obj_classes,
                        target_classes)[:self.num_results]
                    try:
                        npclassscore = full_score[0][0]
                    except IndexError:
                        # no score is found due to missing classes for spectra
                        self._target_no_scores.add(target)  # keep track
                        if not self.filter_missing_scores:
                            results[obj][target] = ObjectLink(
                                obj, target, self, full_score)
                    else:
                        if npclassscore > self.cutoff:
                            results[obj][target] = ObjectLink(
                                obj, target, self, full_score)

        # info about spectra/MFs with missing scoring
        len_missing = len(self._target_no_scores)
        if len_missing > 0:
            filter_msg = 'kept'
            if self.filter_missing_scores:
                filter_msg = 'filtered out'
            logger.warning(
                f'{len_missing} targets have no NPClassScore '
                f'prediction due to missing class predictions and are '
                f'{filter_msg} by default. Adjust .filter_missing_scores '
                f'to change.')

        logger.info(f"NPClassScore completed in {time.time() - begin:.1f}s")
        link_collection._add_links_from_method(self, results)
        return link_collection

    def format_data(self, data):
        """Given whatever output data the method produces, return a readable string version"""
        # data or full_score is a list of tuples, here return just NPClassScore
        formatted_data = None  # default when there is no score (missing class)
        if data:
            # there is a score
            formatted_data = '{:.3f}'.format(data[0][0])
        return formatted_data

    def sort(self, objects, reverse=True):
        """Given a list of objects, return them sorted by link score"""
        return sorted(objects,
                      key=lambda objlink: objlink[self],
                      reverse=reverse)
