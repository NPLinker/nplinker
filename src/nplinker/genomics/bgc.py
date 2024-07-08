from __future__ import annotations
import logging
from typing import TYPE_CHECKING
from deprecated import deprecated
from nplinker.strain import Strain
from .aa_pred import predict_aa


if TYPE_CHECKING:
    from .gcf import GCF

logger = logging.getLogger(__name__)


class BGC:
    """Class to model BGC (biosynthetic gene cluster) data.

    BGC data include both annotations and sequence data. This class is
    mainly designed to model the annotations or metadata.

    The raw BGC data is stored in GenBank format (.gbk). Additional
    [GenBank features](https://www.insdc.org/submitting-standards/feature-table/)
    could be added to the GenBank file to annotate
    BGCs, e.g. antiSMASH has some self-defined features (like `region`) in
    its output GenBank files.

    The annotations of BGC can be stored in JSON format, which is defined
    and used by [MIBiG](https://mibig.secondarymetabolites.org/).

    Attributes:
        id: BGC identifier, e.g. MIBiG accession, GenBank accession.
        product_prediction: A tuple of (predicted) natural
            products or product classes of the BGC.
            For antiSMASH's GenBank data, the feature `region /product`
            gives product information.
            For MIBiG metadata, its biosynthetic class provides such info.
        mibig_bgc_class: A tuple of MIBiG biosynthetic classes to which the BGC belongs.
            Defaults to None.

            MIBiG defines 6 major biosynthetic classes for natural products,
            including `NRP`, `Polyketide`, `RiPP`, `Terpene`, `Saccharide`
            and `Alkaloid`. Note that natural products created by the other
            biosynthetic mechanisms fall under the category `Other`. For more details
            see [the paper](https://doi.org/10.1186/s40793-018-0318-y).
        description: Brief description of the BGC.
            Defaults to None.
        smiles: A tuple of SMILES formulas of the BGC's
            products.
            Defaults to None.
        antismash_file: The path to the antiSMASH GenBank file.
            Defaults to None.
        antismash_id: Identifier of the antiSMASH BGC, referring
            to the feature `VERSION` of GenBank file.
            Defaults to None.
        antismash_region: AntiSMASH BGC region number, referring
            to the feature `region` of GenBank file.
            Defaults to None.
        parents: The set of GCFs that contain the BGC.
        strain: The strain of the BGC.
    """

    def __init__(self, id: str, /, *product_prediction: str):
        """Initialize the BGC object.

        Args:
            id: BGC identifier, e.g. MIBiG accession, GenBank accession.
            product_prediction: BGC's (predicted) natural products or product classes.

        Examples:
            >>> bgc = BGC("Unique_BGC_ID", "Polyketide", "NRP")
            >>> bgc.id
            'Unique_BGC_ID'
            >>> bgc.product_prediction
            ('Polyketide', 'NRP')
            >>> bgc.is_mibig()
            False
        """
        # BGC metadata
        self.id = id
        self.product_prediction = product_prediction

        self.mibig_bgc_class: tuple[str] | None = None
        self.description: str | None = None
        self.smiles: tuple[str] | None = None

        # antismash related attributes
        self.antismash_file: str | None = None
        self.antismash_id: str | None = None  # version in .gbk, id in SeqRecord
        self.antismash_region: int | None = None  # antismash region number

        # other attributes
        self.parents: set[GCF] = set()
        self._strain: Strain | None = None

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "{}(id={}, strain={}, asid={}, region={})".format(
            self.__class__.__name__,
            self.id,
            self.strain,
            self.antismash_id,
            self.antismash_region,
        )

    def __eq__(self, other) -> bool:
        if isinstance(other, BGC):
            return self.id == other.id and self.product_prediction == other.product_prediction
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self.id, self.product_prediction))

    def __reduce__(self) -> tuple:
        """Reduce function for pickling."""
        return (self.__class__, (self.id, *self.product_prediction), self.__dict__)

    def add_parent(self, gcf: GCF) -> None:
        """Add a parent GCF to the BGC.

        Args:
            gcf: gene cluster family
        """
        gcf.add_bgc(self)

    def detach_parent(self, gcf: GCF) -> None:
        """Remove a parent GCF."""
        gcf.detach_bgc(self)

    @property
    def strain(self) -> Strain | None:
        """Get the strain of the BGC."""
        return self._strain

    @strain.setter
    def strain(self, strain: Strain) -> None:
        self._strain = strain

    @property
    def bigscape_classes(self) -> set[str | None]:
        """Get BiG-SCAPE's BGC classes.

        BiG-SCAPE's BGC classes are similar to those defined in MiBIG but have
        more categories (7 classes). More details see:
        https://doi.org/10.1038%2Fs41589-019-0400-9.
        """
        return {p.bigscape_class for p in self.parents}

    def is_mibig(self) -> bool:
        """Check if the BGC is a MIBiG reference BGC or not.

        Warning:
            This method evaluates MIBiG BGC based on the pattern that MIBiG
            BGC names start with "BGC". It might give false positive result.

        Returns:
            True if it's MIBiG reference BGC
        """
        return self.id.startswith("BGC")

    # CG: why not providing whole product but only amino acid as product monomer?
    # this property is not used in NPLinker core business.
    @property
    @deprecated(version="2.0.0", reason="This method will be removed soon")
    def aa_predictions(self) -> list:
        """Amino acids as predicted monomers of product.

        Returns:
            list of dicts with key as amino acid and value as prediction
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
