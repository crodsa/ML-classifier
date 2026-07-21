"""Logging setup and verbose printing helpers."""

import logging


def setup_logging(verbose: bool = False) -> None:
    """Configure root logging level and format.

    Args:
        verbose: If True, sets level to INFO; otherwise WARNING.
    """
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        force=True,
    )


def vprint(msg: str, verbose: bool = False) -> None:
    """Log msg at INFO level when verbose is True.

    Args:
        msg:     Message string.
        verbose: Gate — only logs when True.
    """
    if verbose:
        logging.info(msg)
