from __future__ import annotations
import csv
from os import PathLike
from nplinker.logconfig import LogConfig
from ..abc import GCFLoaderBase
from ..gcf import GCF


logger = LogConfig.getLogger(__name__)


class BigscapeGCFLoader():

    def __init__(self, cluster_file: str | PathLike, /) -> None:
        """Build a loader for BiG-SCAPE GCF cluster file.

        Args:
            cluster_file(str | PathLike): Path to the BiG-SCAPE cluster file,
                the filename has a pattern of "<class>_clustering_c0.xx.tsv".

        Attributes:
            cluster_file(str): path to the BiG-SCAPE clsuter file.
        """
        self.cluster_file = str(cluster_file)
        self._gcf_dict = self._parse_gcf(self.cluster_file)
        self._gcf_list = list(self._gcf_dict.values())

    def get_gcfs(self) -> list[GCF]:
        """Get all GCF objects."""
        return self._gcf_list

    @staticmethod
    def _parse_gcf(cluster_file: str) -> dict[str, GCF]:
        """Parse BiG-SCAPE cluster file to return GCF objects."""
        gcf_dict = {}
        with open(cluster_file, "rt", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # skip headers
            for line in reader:
                bgc_id, family_id = line[:]
                if family_id not in gcf_dict:
                    gcf_dict[family_id] = GCF(family_id)
                gcf_dict[family_id].bgc_ids.add(bgc_id)
        return gcf_dict


# register as virtual class to prevent metaclass conflicts
GCFLoaderBase.register(BigscapeGCFLoader)
