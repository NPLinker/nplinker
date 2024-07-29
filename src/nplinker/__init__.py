import logging
from .logger import setup_logging
from .nplinker import NPLinker


logging.getLogger(__name__).addHandler(logging.NullHandler())

__author__ = "Cunliang Geng"
__email__ = "c.geng@esciencecenter.nl"
__version__ = "2.0.0-alpha.4"


__all__ = ["NPLinker", "setup_logging"]
