import itertools
from nplinker.logconfig import LogConfig


logger = LogConfig.getLogger(__file__)


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
