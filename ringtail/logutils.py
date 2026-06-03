#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail logging module
#

"""
Central logger for the ringtail package.

Usage in any submodule:
    from .logutils import get_logger
    logger = get_logger(__name__)
    logger.debug("message")
    logger.warning("something unexpected: %s", value)

To configure output in an application or script:
    from ringtail import setup_logging
    setup_logging(level="INFO")                      # pretty console output
    setup_logging(level="DEBUG", logfile="rt.log")   # also write to file
"""

import logging
import sys
from typing import Union

LOGGER_NAME = "ringtail"

_root_logger = logging.getLogger(LOGGER_NAME)
_root_logger.addHandler(
    logging.NullHandler()
)  # silent by default (library best practice)


def get_logger(name: str) -> logging.Logger:
    """Return a child logger. Pass __name__ from the calling module."""
    return logging.getLogger(name)


def setup_logging(
    level: Union[str, int] = "WARNING",
    logfile: Union[str, None] = None,
    fmt: str = "%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    datefmt: str = "%H:%M:%S",
) -> None:
    """Configure console (and optionally file) output for the package logger.

    Call once at application startup. Safe to call multiple times —
    existing handlers are replaced each call.
    """
    root = logging.getLogger(LOGGER_NAME)
    root.setLevel(level)
    root.handlers.clear()

    formatter = logging.Formatter(fmt, datefmt=datefmt)

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    root.addHandler(console)

    if logfile:
        fh = logging.FileHandler(logfile, encoding="utf-8")
        fh.setFormatter(formatter)
        root.addHandler(fh)


# Alias so `from ringtail import LOGGER` still works in existing code.
# LOGGER is a standard logging.Logger — use .debug(), .info(), etc. directly.
LOGGER = _root_logger
