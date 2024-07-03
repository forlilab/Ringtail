#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail custom exceptions
#

from .logutils import LOGGER as logger


class StorageError(Exception):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class DatabaseInsertionError(StorageError):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class DatabaseConnectionError(StorageError):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class DatabaseTableCreationError(StorageError):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class DatabaseQueryError(StorageError):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class DatabaseViewCreationError(StorageError):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class RTCoreError(Exception):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class OptionError(Exception):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class FileParsingError(Exception):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class WriteToStorageError(Exception):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class MultiprocessingError(Exception):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class ResultsProcessingError(Exception):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)


class OutputError(Exception):
    def __init__(self, e):
        logger.error(__name__ + ":" + e)
