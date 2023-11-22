from __future__ import annotations
from typing import TYPE_CHECKING
from nplinker.logconfig import LogConfig
from nplinker.strain_collection import StrainCollection


if TYPE_CHECKING:
    from nplinker.strain import Strain
    from .bgc import BGC

logger = LogConfig.getLogger(__name__)


class GCF:
    def __init__(self, gcf_id: str, /) -> None:
        """Class to model gene cluster family (GCF).

        GCF is a group of similar BGCs and generated by clustering BGCs with
        tools such as BiG-SCAPE and BiG-SLICE.

        Args:
            gcf_id(str): id of the GCF object.

        Attributes:
            gcf_id(str): id of the GCF object.
            bgc_ids(set[str]): a set of BGC ids that belongs to the GCF.
            bigscape_class(str | None): BiG-SCAPE's BGC class.
                BiG-SCAPE's BGC classes are similar to those defined in MiBIG
                but have more categories (7 classes). More details see:
                https://doi.org/10.1038%2Fs41589-019-0400-9.
        """
        self.gcf_id = gcf_id
        self.bgc_ids: set[str] = set()
        self.bigscape_class: str | None = None
        self._bgcs: set[BGC] = set()
        self._strains: StrainCollection = StrainCollection()

    def __str__(self) -> str:
        return (
            f"GCF(id={self.gcf_id}, #BGC_objects={len(self.bgcs)}, #bgc_ids={len(self.bgc_ids)},"
            f"#strains={len(self._strains)})."
        )

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: GCF) -> bool:
        if isinstance(other, GCF):
            return self.gcf_id == other.gcf_id and self.bgcs == other.bgcs
        return NotImplemented

    def __hash__(self) -> int:
        """Hash function for GCF.

        Note that GCF class is a mutable container. We only hash the GCF id to
        avoid the hash value changes when `self._bgcs` is updated.
        """
        return hash(self.gcf_id)

    @property
    def bgcs(self) -> set[BGC]:
        """Get the BGC objects."""
        return self._bgcs

    @property
    def strains(self) -> StrainCollection:
        """Get the strains in the GCF."""
        return self._strains

    def add_bgc(self, bgc: BGC) -> None:
        """Add a BGC object to the GCF."""
        bgc.parents.add(self)
        self._bgcs.add(bgc)
        self.bgc_ids.add(bgc.bgc_id)
        if bgc.strain is not None:
            self._strains.add(bgc.strain)
        else:
            logger.warning("No strain specified for the BGC %s", bgc.bgc_id)

    def detach_bgc(self, bgc: BGC) -> None:
        """Remove a child BGC object."""
        bgc.parents.remove(self)
        self._bgcs.remove(bgc)
        self.bgc_ids.remove(bgc.bgc_id)
        if bgc.strain is not None:
            for other_bgc in self._bgcs:
                if other_bgc.strain == bgc.strain:
                    return
            self._strains.remove(bgc.strain)

    def has_strain(self, strain: Strain) -> bool:
        """Check if the given strain exists.

        Args:
            strain(Strain): `Strain` object.

        Returns:
            bool: True when the given strain exist.
        """
        return strain in self._strains

    def has_mibig_only(self) -> bool:
        """Check if the GCF's children are only MIBiG BGCs.

        Returns:
            bool: True if `GCF.bgc_ids` are only MIBiG BGC ids.
        """
        return all(map(lambda id: id.startswith("BGC"), self.bgc_ids))

    def is_singleton(self) -> bool:
        """Check if the GCF contains only one BGC.

        Returns:
            bool: True if `GCF.bgc_ids` contains only one BGC id.
        """
        return len(self.bgc_ids) == 1
