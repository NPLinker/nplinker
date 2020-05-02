import itertools
import random
from types import SimpleNamespace

import numpy as np

from .data_linking import DataLinks, LinkFinder
from ..genomics import BGC, GCF
from ..metabolomics import Spectrum, MolecularFamily

from ..logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class LinkCollection(object):
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
        # object_links is a dict of <source: [ObjectLinks]>
        if method in self._methods:
            # this is probably an error...
            raise Exception('Duplicate method found in LinkCollection: {}'.format(method.name))

        # if this is the first set of results to be generated, can just dump
        # them all straight in
        if len(self._methods) == 0:
            self._link_data = {k: v for k, v in object_links.items()}
        else:
            # if already some results added, in OR mode can just merge the new set
            # with the existing set, but in AND mode need to ensure we end up with
            # only results that appear in both sets
            
            if not self._and_mode:
                self._merge_or_mode(object_links)
            else:
                self._merge_and_mode(object_links)

        self._methods.add(method)

    def _merge_and_mode(self, object_links):
        # set of ObjectLinks common to existing + new results 
        intersect1 = self._link_data.keys() & object_links.keys()

        # iterate over the existing set of link info, remove entries for objects
        # that aren't common to both that and the new set of info, and merge in
        # any common links
        to_remove = []
        for source, existing_links in self._link_data.items():
            if source not in intersect1:
                to_remove.append(source)
                continue

            links_to_merge = object_links[source]
            intersect2 = existing_links.keys() & links_to_merge.keys()
            self._link_data[source] = {k: v for k, v in existing_links.items() if k in intersect2}

            for target, object_link in object_links[source].items():
                self._link_data[source][target]._merge(object_link)

        for source in to_remove:
            del self._link_data[source]

    def _merge_or_mode(self, object_links):
        for source, links in object_links.items():

            # update the existing dict with the new entries that don't appear in it already
            self._link_data[source].update({k: v for k, v in links.items() if k not in self._link_data[source]})

            # now merge the remainder (common to both)
            for target, object_link in links.items():
                self._link_data[source][target]._merge(object_link)

    def filter_no_shared_strains(self):
        len_before = len(self._link_data)
        self.filter_links(lambda x: len(x.shared_strains) > 0)
        logger.debug('filter_no_shared_strains: {} => {}'.format(len_before, len(self._link_data)))

    def filter_sources(self, callable_obj):
        len_before = len(self._link_data)
        self._link_data = {k: v for k, v in self._link_data.items() if callable_obj(k)}
        logger.debug('filter_sources: {} => {}'.format(len_before, len(self._link_data)))

    def filter_targets(self, callable_obj, sources=None):
        to_remove = []
        sources_list = self._link_data.keys() if sources is None else sources
        for source in sources_list:
            self._link_data[source] = {k: v for k, v in self._link_data[source].items() if callable_obj(k)}
            # if there are now no links for this source, remove it completely
            if len(self._link_data[source]) == 0:
                to_remove.append(source)

        for source in to_remove:
            del self._link_data[source]

    def filter_links(self, callable_obj, sources=None):
        to_remove = []
        sources_list = self._link_data.keys() if sources is None else sources
        for source in sources_list:
            self._link_data[source] = {k: v for k, v in self._link_data[source].items() if callable_obj(v)}
            # if there are now no links for this source, remove it completely
            if len(self._link_data[source]) == 0:
                to_remove.append(source)

        for source in to_remove:
            del self._link_data[source]

    def get_sorted_links(self, method, source, reverse=True):
        links = [link for link in self._link_data[source].values() if method in link.methods]
        return method.sort(links, reverse)

    def get_all_targets(self):
        return list(set(itertools.chain.from_iterable(self._link_data[x].keys() for x in self._link_data.keys())))

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

class ObjectLink(object):
    """
    Class which stores information about a single link between two objects.

    The information stored is basically:
     - the "source" of the link (the original object provided as part of the input)
     - the "target" of the link (the linked object)
     - a possibly empty list of Strain objects shared between source and target
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

    @property
    def method_count(self):
        return len(self._method_data)

    @property
    def methods(self):
        return self._method_data.keys()

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
        return 'ObjectLink(source={}, target={}, #methods={})'.format(self.source, self.target, len(self._method_data))

    def __repr__(self):
        return str(self)

class ScoringMethod(object):

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
        results = {obj: data for obj, data in list(mc_results.links.items())[:num_to_keep]}
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

    def get_links(self, objects, link_collection):
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
            type1_objects = objects if spectra_input else {self.npl.spectra[int(index)] for index in results[0][self.R_DST_ID]}

            # other objs should be "objects" if input is GCF, otherwise create 
            # from the results array as above
            other_objs = objects if not spectra_input else {self.npl._gcfs[int(index)] for index in results[0][self.R_DST_ID]}

            # build a lookup table for pairs of object IDs and their index in the results array
            if spectra_input:
                pairs_lookup = {(results[0][self.R_SRC_ID][i], results[0][self.R_DST_ID][i]): i for i in range(len(results[0][self.R_SRC_ID]))}
            else:
                pairs_lookup = {(results[0][self.R_DST_ID][i], results[0][self.R_SRC_ID][i]): i for i in range(len(results[0][self.R_SRC_ID]))}

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
                    final_score = (results[0][self.R_SCORE][k] - expected) / variance_sqrt

                    # finally apply the scoring cutoff and store the result
                    if self.cutoff is None or (final_score >= self.cutoff):
                        new_src.append(results[0][self.R_SRC_ID][k])
                        new_dst.append(results[0][self.R_DST_ID][k])
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
                    if type_ != MolecularFamily:
                        logger.debug('Found no links for {} input objects (type {})'.format(len(objects), type_))
                    continue # no results

                # for each entry in the results (each Spectrum or MolecularFamily)
                for j in range(res.shape[1]):
                    # extract the ID of the object and get the object itself
                    obj_id = int(res[self.R_DST_ID, j])
                    obj = self.npl._spectra[obj_id] if type_ == Spectrum else self.npl._molfams[obj_id]

                    # retrieve the GCF object too (can use its internal ID to index
                    # directly into the .gcfs list)
                    gcf = self.npl._gcfs[int(res[self.R_SRC_ID][j])]

                    # record that this GCF has at least one link associated with it
                    scores_found.add(gcf)

                    # save the scores
                    if gcf not in metcalf_results:
                        metcalf_results[gcf] = {obj: []}
                    metcalf_results[gcf][obj] = ObjectLink(gcf, obj,self, res[self.R_SCORE, j])

        else:
            logger.debug('MetcalfScoring: input_type=Spec/MolFam, result_type=GCF, inputs={}, results={}'.format(len(objects), results[0].shape))
            # for non-GCF input, result is a list containing a single array, shape (3, x)
            # where x is the total number of links found
            results = results[0]
            if results.shape[1] == 0:
                logger.debug('Found no links for {} input objects'.format(len(objects)))
                return link_collection

            # for each entry in the results (each GCF)
            for j in range(results.shape[1]):
                # extract the ID of the GCF and use that to get the object itself
                gcf = self.npl._gcfs[int(results[self.R_DST_ID, j])]

                # retrieve the Spec/MolFam object too (can use its internal ID to index
                # directly into the appropriate list)
                obj_id = int(results[self.R_SRC_ID, j])
                obj = self.npl._spectra[obj_id] if input_type == Spectrum else self.npl._molfams[obj_id]

                # record that this Spectrum or MolecularFamily has at least one link associated with it
                scores_found.add(obj)

                # save the scores
                if obj not in metcalf_results:
                    metcalf_results[obj] = {gcf: []}
                metcalf_results[obj][gcf] = ObjectLink(obj, gcf, self, results[self.R_SCORE, j])

        logger.debug('MetcalfScoring found {} results'.format(len(metcalf_results)))
        link_collection._add_links_from_method(self, metcalf_results)
        logger.debug('MetcalfScoring: completed')
        return link_collection

    def format_data(self, data):
        # for metcalf the data will just be a floating point value (i.e. the score)
        return '{:.4f}'.format(data)

    def sort(self, objects, reverse=True):
        # sort based on score
        return sorted(objects, key=lambda objlink: objlink[self], reverse=reverse)
