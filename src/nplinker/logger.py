import logging
from rich.console import Console
from rich.logging import RichHandler


def setup_logging(level: str = "INFO", file: str = "", use_console: bool = True) -> None:
    """Setup logging configuration for the ancestor logger "nplinker".

    ??? info "Usage Documentation"
        [How to setup logging][how-to-setup-logging]

    Args:
        level: The log level, use the logging module's log level constants.
            Valid levels are: `NOTSET`, `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`.
        file: The file to write the log to.
            If the file is an empty string (by default), the log will not be written to a file.
            If the file does not exist, it will be created.
            The log will be written to the file in append mode.
        use_console: Whether to log to the console.
    """
    # Get the ancestor logger "nplinker"
    logger = logging.getLogger("nplinker")
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
