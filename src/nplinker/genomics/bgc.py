from __future__ import annotations
import logging
from typing import TYPE_CHECKING
from typing import Any
from nplinker.strain import Strain


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
            Defaults to None, which means the class is unknown.

            MIBiG defines 6 major biosynthetic classes for natural products,
            including `NRP`, `Polyketide`, `RiPP`, `Terpene`, `Saccharide`
            and `Alkaloid`. Note that natural products created by the other
            biosynthetic mechanisms fall under the category `Other`. For more details
            see [the paper](https://doi.org/10.1186/s40793-018-0318-y).
        description: Brief description of the BGC.
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
        more categories (7 classes), including:

        - NRPS
        - PKS-NRP_Hybrids
        - PKSI
        - PKSother
        - RiPPs
        - Saccharides
        - Terpene

        For BGC falls outside of these categories, the value is "Others".

        Default is None, which means the class is unknown.

        More details see:
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

    def to_dict(self) -> dict[str, Any]:
        """Convert the BGC object to a dictionary for exporting purpose.

        Returns:
            A dictionary containing the following key-value pairs:

            - GCF_id (list[str]): A list of GCF IDs.
            - GCF_bigscape_class (list[str]): A list of BiG-SCAPE classes.
            - strain_id (str | None): The ID of the strain.
            - description (str | None): A description of the BGC.
            - BGC_name (str): The name of the BGC.
            - product_prediction (list[str]): (predicted) products or product classes of the BGC.
            - mibig_bgc_class (list[str] | None): MIBiG biosynthetic classes.
            - antismash_id (str | None): The antiSMASH ID.
            - antismash_region (int | None): The antiSMASH region number.
        """
        # Keys are ordered to make the output easier to analyze
        return {
            "GCF_id": [gcf.id for gcf in self.parents if gcf.id is not None],
            "GCF_bigscape_class": [bsc for bsc in self.bigscape_classes if bsc is not None],
            "strain_id": self.strain.id if self.strain is not None else None,
            "description": self.description,
            "BGC_name": self.id,
            "product_prediction": list(self.product_prediction),
            "mibig_bgc_class": self.mibig_bgc_class,
            "antismash_id": self.antismash_id,
            "antismash_region": self.antismash_region,
        }

    def to_tabular(self) -> dict[str, str]:
        """Convert the BGC object to a tabular format.

        Returns:
            dict: A dictionary representing the BGC object in tabular format.
                The keys can be treated as headers and values are strings in which tabs are removed.
                This dict can be exported as a TSV file.
        """
        return {
            key: self._to_string(value).replace("\t", "    ")
            for key, value in self.to_dict().items()
        }

    @staticmethod
    def _to_string(value: Any) -> str:
        """Convert various types of values to a string.

        Args:
            value: The value to be converted to a string.
                Can be a list, dict, or any other JSON-compatible type.

        Returns:
            A string representation of the input value.
        """
        # Convert list to comma-separated string
        if isinstance(value, list):
            formatted_value = ", ".join(map(str, value))
        # Convert dict to comma-separated string
        elif isinstance(value, dict):
            formatted_value = ", ".join([f"{k}:{v}" for k, v in value.items()])
        # Convert None to empty string
        elif value is None:
            formatted_value = ""
        # Convert anything else to string
        else:
            formatted_value = str(value)
        return formatted_value
