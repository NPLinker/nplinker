from __future__ import annotations
import fnmatch
import os
from Bio import SeqIO
from Bio import SeqRecord
from nplinker.genomics import BGC
from nplinker.logconfig import LogConfig
from nplinker.strains import Strain
from nplinker.utils import list_dirs
from nplinker.utils import list_files
from ..abc import BGCLoaderBase


logger = LogConfig.getLogger(__name__)


class AntismashBGCLoader:

    def __init__(self, data_dir: str) -> None:
        """Build a loader for AntiSMASH BGC genbank (.gbk) files

        Note:
            AntiSMASH BGC directory must follow the structure below:
            antismash
                ├── antismash_id_1 (one AntiSMASH output, e.g. GCF_000514775.1)
                │  ├── GCF_000514775.1.gbk
                │  ├── NZ_AZWO01000004.region001.gbk
                │  └── ...
                ├── antismash_id_2
                │  ├── ...
                └── ...

        Args:
            antismash_dir(str): Path to AntiSMASH directory that contains a
                collection of AntiSMASH outputs.
        """
        self.data_dir = data_dir
        self._file_dict = self._parse_data_dir(self.data_dir)
        self._bgc_dict = self._parse_bgcs(self._file_dict)

    def get_files(self) -> dict[str, str]:
        """Get BGC gbk files

        Returns:
            dict[str, str]: key is BGC name (gbk file name) and value is path to
                the gbk file
        """
        return self._file_dict

    @staticmethod
    def _parse_data_dir(data_dir: str) -> dict[str, str]:
        """Parse AntiSMASH directory to get path of all BGC gbk files.

        Args:
            data_dir(str): Path to AntiSMASH directory that contains
                a collection of AntiSMASH outputs

        Returns:
            dict[str, str]: key is BGC name (gbk file name) and value is path to
                the gbk file
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

    def get_bgcs(self) -> dict[str, BGC]:
        """Get all BGC objects

        Returns:
            dict[str, BGC]: key is BGC name and value is
                :class:`~nplinker.genomic.BGC` objects
        """
        return self._bgc_dict

    @staticmethod
    def _parse_bgcs(bgc_files: dict[str, str]) -> dict[str, BGC]:
        """Load given BGC files as BGC objects

        Args:
            bgc_files(dict[str, str]): key is BGC name and value is path to the
                BGC gbk file, see method :meth:`.bgc_files`.

        Returns:
            dict[str, BGC]: key is BGC name and value is :class:`~nplinker.genomic.BGC` objects
        """
        bgcs = {}
        for bgc_id in bgc_files:
            bgc = parse_bgc_genbank(bgc_files[bgc_id])
            bgcs[bgc_id] = bgc
        return bgcs


def parse_bgc_genbank(file: str) -> BGC:
    """Parse a single BGC gbk file to BGC object.

    Note:
        If product info is not available in gbk file, the product of BGC
            object (bgc.product_prediction) is set to empty list.

    Args:
        file(str): Path to BGC gbk file

    Returns:
        BGC: :class:`~nplinker.genomic.BGC` object

    Examples:
        >>> bgc = AntismashBGCLoader.parse_bgc(
        ...    "/data/antismash/GCF_000016425.1/NC_009380.1.region001.gbk")
    """
    fname = os.path.splitext(os.path.basename(file))[0]

    record = SeqIO.read(file, format="genbank")
    description = record.description  # "DEFINITION" in gbk file
    antismash_id = record.id  # "VERSION" in gbk file
    features = _parse_antismash_genbank(record)
    product_prediction = features.get("product")
    if product_prediction is None:
        raise ValueError(
            f"Not found product prediction in antiSMASH Genbank file {file}")

    # init BGC
    bgc = BGC(bgc_id=fname, product_prediction=product_prediction)
    bgc.description = description
    bgc.antismash_id = antismash_id
    bgc.antismash_file = file
    bgc.antismash_region = features.get("region_number")
    bgc.smiles = features.get("smiles")
    bgc.strain = Strain(fname)
    return bgc


def _parse_antismash_genbank(record: SeqRecord.SeqRecord) -> dict:
    features = {}
    for feature in record.features:
        if feature.type == "region":
            # biopython assumes region numer is a list, but it's actually an int
            features["region_number"] = feature.qualifiers.get('region_number')[0]
            features["product"] = feature.qualifiers.get('product')
        if feature.type == "cand_cluster":
            smiles = feature.qualifiers.get('SMILES')
            # space is not allowed in SMILES spec
            # biopython generates space when reading multi-line SMILES from .gbk
            if smiles is not None:
                smiles = [i.replace(' ', '') for i in smiles]
            features["smiles"] = smiles
    return features


# register as virtual class to prevent metaclass conflicts
BGCLoaderBase.register(AntismashBGCLoader)
