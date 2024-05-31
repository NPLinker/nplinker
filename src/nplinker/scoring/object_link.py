class ObjectLink:
    """Class which stores information about a single link between two objects.

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

    def __init__(self, source, target, method, data=None, common_strains=[]):
        self.source = source
        self.target = target
        self.common_strains = common_strains
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
        # TODO: hashable object should also have `__eq__` defined, see #136.
        # this implementation is not ideal as the hash value is not unique
        return hash(self.source.id)

    def __str__(self):
        return "ObjectLink(source={}, target={}, #methods={})".format(
            self.source, self.target, len(self._method_data)
        )

    def __repr__(self):
        return str(self)
