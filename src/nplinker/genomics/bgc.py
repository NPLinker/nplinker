from __future__ import annotations
from typing import TYPE_CHECKING
from deprecated import deprecated
from nplinker.logconfig import LogConfig
from .aa_pred import predict_aa

if TYPE_CHECKING:
    from ..strains import Strain
    from .gcf import GCF

logger = LogConfig.getLogger(__name__)


class BGC():

    def __init__(self,
                 bgc_id: str,
                 product_prediction: list[str]):
        self.bgc_id = bgc_id  # BGC file name
        self.product_prediction = product_prediction  # can get from gbk SeqFeature "region"

        # MIBiG biosynthetic class
        self.mibig_bgc_class: list[str] | None = None

        self.description: str | None = None

        # CG TODO: change parents to parent
        self.parents: set[GCF] = set()
        self.smiles: list[str] = []

        # antismash related
        self.antismash_file: str | None = None
        self.antismash_id: str | None = None  # version in .gbk, id in SeqRecord
        self.antismash_region: int | None = None  # antismash region number

        self._strain: Strain | None = None

    def add_parent(self, gcf):
        self.parents.add(gcf)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return '{}(bgc_id={}, strain={}, asid={}, region={})'.format(
            self.__class__.__name__, self.bgc_id, self.strain,
            self.antismash_id, self.antismash_region)

    def __eq__(self, other):
        return self.bgc_id == other.bgc_id

    def __hash__(self):
        return self.bgc_id

    @property
    def strain(self) -> Strain | None:
        return self._strain

    @strain.setter
    def strain(self, strain: Strain) -> None:
        self._strain = strain

    @property
    def bigscape_classes(self):
        return {p.bigscape_class for p in self.parents}

    def is_mibig(self):
        """Check if the BGC is MIBiG reference BGC or not.

        Note:
            This method evaluates MIBiG BGC based on the pattern that MIBiG
            BGC names start with "BGC". It might give false positive result.

        Returns:
            bool: True if it's MIBiG reference BGC
        """
        return self.bgc_id.startswith('BGC')

    # CG: why not providing whole product but only amino acid as product monomer?
    # this property is not used in NPLinker core business.
    @property
    @deprecated(version='2.0.0', reason="This method will be removed soon")
    def aa_predictions(self):
        """Amino acids as predicted monomers of product.

        Returns:
            list: list of dicts with key as amino acid and value as prediction
                probability.
        """
        # Load aa predictions and cache them
        self._aa_predictions = None
        if self._aa_predictions is None:
            self._aa_predictions = {}
            if self.antismash_file is not None:
                for p in predict_aa(self.antismash_file):
                    self._aa_predictions[p[0]] = p[1]
        return [self._aa_predictions]
