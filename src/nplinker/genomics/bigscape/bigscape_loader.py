from __future__ import annotations
import csv
import logging
import sqlite3
from os import PathLike
from ..abc import GCFLoaderBase
from ..gcf import GCF


logger = logging.getLogger(__name__)


class BigscapeGCFLoader(GCFLoaderBase):
    """Data loader for BiG-SCAPE GCF cluster file.

    Attributes:
        cluster_file: path to the BiG-SCAPE cluster file.
    """

    def __init__(self, cluster_file: str | PathLike, /) -> None:
        """Initialize the BiG-SCAPE GCF loader.

        Args:
            cluster_file: Path to the BiG-SCAPE cluster file,
                the filename has a pattern of `<class>_clustering_c0.xx.tsv`.
        """
        self.cluster_file: str = str(cluster_file)
        self._gcf_list = self._parse_gcf(self.cluster_file)

    def get_gcfs(self, keep_mibig_only: bool = False, keep_singleton: bool = False) -> list[GCF]:
        """Get all GCF objects.

        Args:
            keep_mibig_only: True to keep GCFs that contain only MIBiG
                BGCs.
            keep_singleton: True to keep singleton GCFs. A singleton GCF
                is a GCF that contains only one BGC.

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


class BigscapeV2GCFLoader(GCFLoaderBase):
    """Data loader for BiG-SCAPE v2 database file.

    Attributes:
        db_file: Path to the BiG-SCAPE database file.
    """

    def __init__(self, db_file: str | PathLike, /) -> None:
        """Initialize the BiG-SCAPE v2 GCF loader.

        Args:
            db_file: Path to the BiG-SCAPE v2 database file
        """
        self.db_file = str(db_file)
        self._gcf_list = self._parse_gcf(self.db_file)

    def get_gcfs(self, keep_mibig_only: bool = False, keep_singleton: bool = False) -> list[GCF]:
        """Get all GCF objects.

        Args:
            keep_mibig_only: True to keep GCFs that contain only MIBiG BGCs.
            keep_singleton: True to keep singleton GCFs.
                A singleton GCF is a GCF that contains only one BGC.

        Returns:
            a list of GCF objects.
        """
        gcf_list = self._gcf_list
        if not keep_mibig_only:
            gcf_list = [gcf for gcf in gcf_list if not gcf.has_mibig_only()]
        if not keep_singleton:
            gcf_list = [gcf for gcf in gcf_list if not gcf.is_singleton()]
        return gcf_list

    @staticmethod
    def _parse_gcf(db_file: str) -> list[GCF]:
        """Get GCF objects from database.

        Args:
            db_file: Path to the sqlite3 database file.

        Returns:
            A list of GCF objects
        """
        gcf_dict: dict[str, GCF] = {}

        with sqlite3.connect(db_file) as connection:
            cursor = connection.cursor()

            query = """
            SELECT gbk.path, bgc_record_family.family_id FROM bgc_record_family
            JOIN bgc_record ON bgc_record.id = bgc_record_family.record_id
            JOIN gbk ON gbk.id = bgc_record.gbk_id
            """

            results = cursor.execute(query).fetchall()

            for result in results:
                gbk_path, family_id = result

                # take the filename of the gbk path as the bgc_id
                # filename
                bgc_id: str = gbk_path.split("/")[-1]
                # remove extension
                bgc_id = bgc_id.rsplit(".", 1)[0]

                if family_id not in gcf_dict:
                    gcf_dict[family_id] = GCF(family_id)
                gcf_dict[family_id].bgc_ids.add(bgc_id)

        return list(gcf_dict.values())
