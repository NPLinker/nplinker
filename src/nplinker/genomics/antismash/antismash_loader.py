from __future__ import annotations
import fnmatch
import logging
import os
from collections.abc import Mapping
from os import PathLike
from pathlib import Path
from Bio import SeqIO
from Bio import SeqRecord
from nplinker.genomics import BGC
from nplinker.strain import Strain
from nplinker.utils import list_dirs
from nplinker.utils import list_files
from ..abc import BGCLoaderBase


logger = logging.getLogger(__name__)


class AntismashBGCLoader(BGCLoaderBase):
    """Data loader for AntiSMASH BGC genbank (.gbk) files."""

    def __init__(self, data_dir: str | PathLike) -> None:
        """Initialize the AntiSMASH BGC loader.

        Args:
            data_dir: Path to AntiSMASH directory that contains a collection of AntiSMASH outputs.

        Notes:
            The input `data_dir` must follow the structure defined in the
            [Working Directory Structure][working-directory-structure] for AntiSMASH data, e.g.:
            ```shell
            antismash
                ├── genome_id_1                  # one AntiSMASH output, e.g. GCF_000514775.1
                │  ├── NZ_AZWO01000004.region001.gbk
                │  └── ...
                ├── genome_id_2
                │  ├── ...
                └── ...
            ```
        """
        self.data_dir = str(data_dir)
        self._file_dict = self._parse_data_dir(self.data_dir)
        self._bgcs = self._parse_bgcs(self._file_dict)

    def get_bgc_genome_mapping(self) -> dict[str, str]:
        """Get the mapping from BGC to genome.

        !!! info
            The directory name of the gbk files is treated as genome id.

        Returns:
            The key is BGC name (gbk file name) and value is genome id (the directory name of the
            gbk file).
        """
        return {
            bid: os.path.basename(os.path.dirname(bpath)) for bid, bpath in self._file_dict.items()
        }

    def get_files(self) -> dict[str, str]:
        """Get BGC gbk files.

        Returns:
            The key is BGC name (gbk file name) and value is path to the gbk file.
        """
        return self._file_dict

    @staticmethod
    def _parse_data_dir(data_dir: str) -> dict[str, str]:
        """Parse AntiSMASH directory to get path of all BGC gbk files.

        Args:
            data_dir: Path to AntiSMASH directory that contains
                a collection of AntiSMASH outputs

        Returns:
            The key is BGC name (gbk file name) and value is path to the gbk file.
        """
        bgc_files = {}
        subdirs = list_dirs(data_dir)
        for subdir in subdirs:
            # get all .gbk files
            files = list_files(subdir, suffix=".gbk", keep_parent=False)
            # filter BGC's .gbk files
            files = fnmatch.filter(files, "*.region???.gbk")
            for f in files:
                fname = os.path.splitext(f)[0]
                fpath = os.path.join(subdir, f)
                bgc_files[fname] = fpath

        return bgc_files

    def get_bgcs(self) -> list[BGC]:
        """Get all BGC objects.

        Returns:
            A list of BGC objects
        """
        return self._bgcs

    @staticmethod
    def _parse_bgcs(bgc_files: Mapping[str, str]) -> list[BGC]:
        """Load given BGC files as BGC objects.

        Args:
            bgc_files: key is BGC name and value is path to the
                BGC gbk file, see method :meth:`.bgc_files`.

        Returns:
            A list of BGC objects
        """
        return [parse_bgc_genbank(file) for file in bgc_files.values()]


def parse_bgc_genbank(file: str | PathLike) -> BGC:
    """Parse a single BGC gbk file to BGC object.

    Args:
        file: Path to BGC gbk file

    Returns:
        BGC object

    Examples:
        >>> bgc = AntismashBGCLoader.parse_bgc(
        ...    "/data/antismash/GCF_000016425.1/NC_009380.1.region001.gbk")
    """
    file = Path(file)
    fname = file.stem

    record = SeqIO.read(file, format="genbank")
    description = record.description  # "DEFINITION" in gbk file
    antismash_id = record.id  # "VERSION" in gbk file
    features = _parse_antismash_genbank(record)
    product_prediction = features.get("product")
    if product_prediction is None:
        raise ValueError(f"Not found product prediction in antiSMASH Genbank file {file}")

    # init BGC
    bgc = BGC(fname, *product_prediction)
    bgc.description = description
    bgc.antismash_id = antismash_id
    bgc.antismash_file = str(file)
    bgc.antismash_region = features.get("region_number")
    bgc.smiles = features.get("smiles")
    bgc.strain = Strain(fname)
    return bgc


def _parse_antismash_genbank(record: SeqRecord.SeqRecord) -> dict:
    features = {}
    for feature in record.features:
        if feature.type == "region":
            # biopython assumes region numer is a list, but it's actually an int
            features["region_number"] = feature.qualifiers.get("region_number")[0]
            features["product"] = feature.qualifiers.get("product")
        if feature.type == "cand_cluster":
            smiles = feature.qualifiers.get("SMILES")
            # space is not allowed in SMILES spec
            # biopython generates space when reading multi-line SMILES from .gbk
            if smiles is not None:
                smiles = tuple(i.replace(" ", "") for i in smiles)
            features["smiles"] = smiles
    return features
