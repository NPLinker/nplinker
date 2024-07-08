from __future__ import annotations
import json
from os import PathLike


class MibigMetadata:
    """Class to model the BGC metadata/annotations defined in MIBiG.

    MIBiG is a specification of BGC metadata and use JSON schema to
    represent BGC metadata. More details see:
    https://mibig.secondarymetabolites.org/download.
    """

    def __init__(self, file: str | PathLike) -> None:
        """Initialize the MIBiG metadata object.

        Args:
            file: Path to the json file of MIBiG BGC metadata

        Examples:
            >>> metadata = MibigMetadata("/data/BGC0000001.json")
        """
        self.file = str(file)
        with open(self.file, "rb") as f:
            self.metadata = json.load(f)

        self._mibig_accession: str
        self._biosyn_class: tuple[str]
        self._parse_metadata()

    @property
    def mibig_accession(self) -> str:
        """Get the value of metadata item 'mibig_accession'."""
        return self._mibig_accession

    @property
    def biosyn_class(self) -> tuple[str]:
        """Get the value of metadata item 'biosyn_class'.

        The 'biosyn_class' is biosynthetic class(es), namely the type of
        natural product or secondary metabolite.

        MIBiG defines 6 major biosynthetic classes for natural products,
        including `NRP`, `Polyketide`, `RiPP`, `Terpene`, `Saccharide`
        and `Alkaloid`. Note that natural products created by the other
        biosynthetic mechanisms fall under the category `Other`. For more details
        see [the paper](https://doi.org/10.1186/s40793-018-0318-y).
        """
        return self._biosyn_class

    def _parse_metadata(self) -> None:
        """Parse metadata to get 'mibig_accession' and 'biosyn_class' values."""
        if "general_params" in self.metadata:
            self._mibig_accession = self.metadata["general_params"]["mibig_accession"]
            self._biosyn_class = tuple(self.metadata["general_params"]["biosyn_class"])
        else:  # version≥2.0
            self._mibig_accession = self.metadata["cluster"]["mibig_accession"]
            self._biosyn_class = tuple(self.metadata["cluster"]["biosyn_class"])
