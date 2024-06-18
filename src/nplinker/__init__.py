import logging
from pathlib import Path
from .logger import setup_logging
from .nplinker import NPLinker


logging.getLogger(__name__).addHandler(logging.NullHandler())

__author__ = "Cunliang Geng"
__email__ = "c.geng@esciencecenter.nl"
__version__ = "2.0.0-alpha.1"


# The path to the NPLinker application database directory
NPLINKER_APP_DATA_DIR = Path(__file__).parent / "data"
del Path


__all__ = ["NPLinker", "setup_logging", "NPLINKER_APP_DATA_DIR"]
