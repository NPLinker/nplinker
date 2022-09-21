import logging


logging.getLogger(__name__).addHandler(logging.NullHandler())

__author__ = "Helge Hecht"
__email__ = "h.hecht@esciencecenter.nl"
__version__ = "1.3.2"

from .BGC import BGC
from .GCF import GCF
from .genomics import loadBGC_from_cluster_files
from .genomics import make_mibig_bgc_dict
from .MiBIGBGC import MiBIGBGC


__all__ = [
    "BGC",
    "MiBIGBGC",
    "GCF",
    "loadBGC_from_cluster_files",
    "make_mibig_bgc_dict"
]
