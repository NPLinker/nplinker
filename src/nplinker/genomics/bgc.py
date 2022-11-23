from __future__ import annotations
import re
from typing import TYPE_CHECKING
from nplinker.logconfig import LogConfig
from .aa_pred import predict_aa
from .genomics_utilities import get_smiles

if TYPE_CHECKING:
    from ..strains import Strain
    from .gcf import GCF

logger = LogConfig.getLogger(__name__)

CLUSTER_REGION_REGEX = re.compile('(.+?)\\.(cluster|region)(\\d+).gbk$')


class BGC():

    def __init__(self,
                 id: int,
                 strain: Strain,
                 name: str,
                 product_prediction: list[str],
                 description: str | None = None):
        self.id = id
        self.strain = strain
        self.name = name  # BGC file name
        self.product_prediction = product_prediction  # can get from gbk SeqFeature "region"
        # allow for multiple parents in the case of hybrid BGCs
        self.parents: set[GCF] = set()
        self.description = description  # can get from gbk SeqRecord.description
        # these will get parsed from the .gbk file
        self.antismash_id: str | None = None  # version in .gbk, id in SeqRecord
        self.antismash_accession: str | None = None  # accession in .gbk, name in SeqRecord

        self.region = -1
        self.cluster = -1

        self.antismash_file = None
        self._aa_predictions = None
        self._known_cluster_blast = None
        self._smiles = None
        self._smiles_parsed = False

        self.edges: set = set()

    def set_filename(self, filename):
        self.antismash_file = filename
        if filename is None:
            return

        # try to extract the region or cluster number if provided
        regex_obj = CLUSTER_REGION_REGEX.search(filename)
        if regex_obj is not None:
            c_or_r = regex_obj.group(2)
            if c_or_r == 'region':
                self.region = int(regex_obj.group(3))
            elif c_or_r == 'cluster':
                self.cluster = int(regex_obj.group(3))

    def add_parent(self, gcf):
        self.parents.add(gcf)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return '{}(id={}, name={}, strain={}, asid={}, region={})'.format(
            self.__class__.__name__, self.id, self.name, self.strain,
            self.antismash_id, self.region)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

    @property
    def bigscape_classes(self):
        return {p.bigscape_class for p in self.parents}

    @property
    def is_hybrid(self):
        return len(self.parents) > 1

    @property
    def is_mibig(self):
        """Check if the BGC is MIBiG reference BGC or not.

        Note:
            This method evaluates MIBiG BGC based on the pattern that MIBiG
            BGC names start with "BGC". It might give false positive result.

        Returns:
            bool: True if it's MIBiG reference BGC
        """
        return self.name.startswith('BGC')

    @property
    def smiles(self):
        if self._smiles is not None or self._smiles_parsed:
            return self._smiles

        if self.antismash_file is None:
            return None

        self._smiles = get_smiles(self)
        self._smiles_parsed = True
        logger.debug(f'SMILES for {self} = {self._smiles}')
        return self._smiles

    # CG: why not providing whole product but only amino acid as product monomer?
    # this property is not used in NPLinker core business.
    @property
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
