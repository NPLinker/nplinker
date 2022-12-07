from nplinker.strain_collection import StrainCollection
from .bgc import BGC


class GCF():
    # CG: the most importance info is family id,
    # the other annotations can be find in BGC gbk files

    def __init__(self, id, gcf_id, bigscape_class):
        self.id = id
        self.gcf_id = gcf_id
        self.bigscape_class = bigscape_class
        self.bgcs = set()

        self.strains = StrainCollection()
        self.strains_lookup = {}

    def __str__(self):
        return 'GCF(id={}, class={}, gcf_id={}, strains={})'.format(
            self.id, self.bigscape_class, self.gcf_id, len(self.strains))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

    def add_bgc(self, bgc):
        self.bgcs.add(bgc)
        bgc.add_parent(self)
        self.strains.add(bgc.strain)
        self.strains_lookup[bgc.strain] = bgc

    def has_strain(self, strain):
        return strain in self.strains
