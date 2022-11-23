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

        self._aa_predictions = None
        self.strains = StrainCollection()
        self.strains_lookup = {}

        # will be set to False if the GCF ends up containing any "hybrid" BGCs
        self._is_pure = True

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
        if bgc.is_hybrid:
            self._is_pure = False

    def has_strain(self, strain):
        return strain in self.strains

    def bgc_for_strain(self, strain):
        return self.strains_lookup[strain]

    def only_mibig(self):
        return len(self.bgcs) == self.num_mibig_bgcs

    def has_mibig(self):
        return self.num_mibig_bgcs > 0

    @property
    def is_pure(self):
        return self._is_pure

    @property
    def non_mibig_bgcs(self):
        return list(
            filter(lambda bgc: not bgc.is_mibig(), self.bgcs))

    @property
    def mibig_bgcs(self):
        return list(filter(lambda bgc: bgc.is_mibig(), self.bgcs))

    @property
    def num_mibig_bgcs(self):
        return len(self.mibig_bgcs)

    @property
    def num_non_mibig_bgcs(self):
        return len(self.bgcs) - self.num_mibig_bgcs

    @property
    def aa_predictions(self):
        """
        Return the predicted AAs for the GCF
        """
        if self._aa_predictions is None:
            bgc_aa_prob = []
            for bgc_count, bgc in enumerate(self.bgcs):
                if not bgc.name.startswith('BGC'):
                    bgc_aa_prob.extend(bgc.aa_predictions)
            self._aa_predictions = bgc_aa_prob

        return self._aa_predictions
