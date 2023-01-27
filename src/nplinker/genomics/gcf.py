from nplinker.strain_collection import StrainCollection
from .bgc import BGC


class GCF():
    def __init__(self, gcf_id: str, bigscape_class: str):
        """Class to represent gene cluster family.

        Args:
            gcf_id(str): family id
            bigscape_class(str): class predicted by bigscape
        """
        self.gcf_id = gcf_id
        self.bigscape_class = bigscape_class
        self.bgcs = set()

        self.strains = StrainCollection()
        self.strains_lookup = {}

    def __str__(self):
        return 'GCF(gcf_id={}, class={},  strains={})'.format(
            self.gcf_id, self.bigscape_class, len(self.strains))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.gcf_id == other.gcf_id

    def __hash__(self):
        return self.gcf_id

    def add_bgc(self, bgc):
        self.bgcs.add(bgc)
        bgc.add_parent(self)
        self.strains.add(bgc.strain)
        self.strains_lookup[bgc.strain] = bgc

    def has_strain(self, strain):
        return strain in self.strains
