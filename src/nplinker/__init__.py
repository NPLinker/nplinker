import logging
from pathlib import Path


logging.getLogger(__name__).addHandler(logging.NullHandler())

__author__ = "Cunliang Geng"
__email__ = "c.geng@esciencecenter.nl"
__version__ = "2.0.0-alpha.1"


# The path to the NPLinker application database directory
NPLINKER_APP_DATA_DIR = Path(__file__).parent / "data"
del Path


def setup_logging(level: str = "INFO", file: str = "", use_console: bool = True) -> None:
    """Setup logging configuration for the ancestor logger "nplinker".

    Args:
        level: The log level, use the logging module's log level constants. Valid levels are:
            "NOTSET", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL".
        file: The file to write the log to. If the file does not exist, it will be created. The log
            will be written to the file in append mode. If the file is an empty string (by default),
            the log will not be written to a file.
        use_console: Whether to log to the console.
    """
    from rich.console import Console
    from rich.logging import RichHandler

    # Get the ancestor logger "nplinker"
    logger = logging.getLogger(__name__)
    logger.setLevel(level)

    # File handler
    if file:
        logger.addHandler(
            RichHandler(
                console=Console(file=open(file, "a"), width=120),  # force the line width to 120
                omit_repeated_times=False,
                rich_tracebacks=True,
                tracebacks_show_locals=True,
                log_time_format="[%Y-%m-%d %X]",
            )
        )

    # Console handler
    if use_console:
        logger.addHandler(
            RichHandler(
                omit_repeated_times=False,
                rich_tracebacks=True,
                tracebacks_show_locals=True,
                log_time_format="[%Y-%m-%d %X]",
            )
        )
