from __future__ import annotations
import csv
from os import PathLike
from nplinker.logconfig import LogConfig
from ..abc import GCFLoaderBase
from ..gcf import GCF


logger = LogConfig.getLogger(__name__)


class BigscapeGCFLoader:
    def __init__(self, cluster_file: str | PathLike, /) -> None:
        """Build a loader for BiG-SCAPE GCF cluster file.

        Args:
            cluster_file: Path to the BiG-SCAPE cluster file,
                the filename has a pattern of "<class>_clustering_c0.xx.tsv".

        Attributes:
            cluster_file: path to the BiG-SCAPE cluster file.
        """
        self.cluster_file = str(cluster_file)
        self._gcf_list = self._parse_gcf(self.cluster_file)

    def get_gcfs(self, keep_mibig_only=False, keep_singleton=False) -> list[GCF]:
        """Get all GCF objects.

        Args:
            keep_mibig_only: True to keep GCFs that contain only MIBiG
                BGCs.
            keep_singleton: True to keep singleton GCFs. A singleton GCF
                is a GCF that contains only one BGC.

        Returns:
            list[GCF]: a list of GCF objects.
        """
        gcf_list = self._gcf_list
        if not keep_mibig_only:
            gcf_list = [gcf for gcf in gcf_list if not gcf.has_mibig_only()]
        if not keep_singleton:
            gcf_list = [gcf for gcf in gcf_list if not gcf.is_singleton()]
        return gcf_list

    @staticmethod
    def _parse_gcf(cluster_file: str) -> list[GCF]:
        """Parse BiG-SCAPE cluster file to return GCF objects."""
        gcf_dict: dict[str, GCF] = {}
        with open(cluster_file, "rt", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # skip headers
            for line in reader:
                bgc_id, family_id = line[:]
                if family_id not in gcf_dict:
                    gcf_dict[family_id] = GCF(family_id)
                gcf_dict[family_id].bgc_ids.add(bgc_id)
        return list(gcf_dict.values())


# register as virtual class to prevent metaclass conflicts
GCFLoaderBase.register(BigscapeGCFLoader)
