from __future__ import annotations
import csv
import logging
from os import PathLike
from .abc import GCFLoaderBase
from .gcf import GCF


logger = logging.getLogger(__name__)


class FakeGCFLoader(GCFLoaderBase):
    """Data loader for Fake GCF cluster file.

    Attributes:
        cluster_file: path to the Fake cluster file.
    """

    def __init__(self, cluster_file: str | PathLike, /) -> None:
        """Initialize the Fake GCF loader.

        Args:
            cluster_file: Path to the Fake cluster file.
        """
        self.cluster_file: str = str(cluster_file)
        self._gcf_list = self._parse_gcf(self.cluster_file)

    def get_gcfs(self, keep_mibig_only: bool = False, keep_singleton: bool = False) -> list[GCF]:
        """Get all GCF objects.

        Args:
            keep_mibig_only: True to keep GCFs that contain only MIBiG BGCs.
            keep_singleton: True to keep singleton GCFs, which are GCFs that contains only one BGC.

        Returns:
            A list of GCF objects.
        """
        gcf_list = self._gcf_list
        if not keep_mibig_only:
            gcf_list = [gcf for gcf in gcf_list if not gcf.has_mibig_only()]
        if not keep_singleton:
            gcf_list = [gcf for gcf in gcf_list if not gcf.is_singleton()]
        return gcf_list

    @staticmethod
    def _parse_gcf(cluster_file: str) -> list[GCF]:
        """Parse Fake cluster file to return GCF objects."""
        gcf_dict: dict[str, GCF] = {}
        with open(cluster_file, "rt", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=",")
            next(reader)  # skip headers
            for line in reader:
                gcf_id, bgc_id = line[:]
                if gcf_id not in gcf_dict:
                    gcf_dict[gcf_id] = GCF(gcf_id)
                gcf_dict[gcf_id].bgc_ids.add(bgc_id)
        return list(gcf_dict.values())
