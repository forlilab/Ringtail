#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail custom exceptions
#

from .logmanager import logger

class StorageError(Exception):
    def __init__(self, e):
        logger.error(e)


class DatabaseInsertionError(StorageError):
    def __init__(self, e):
        logger.error(e)


class DatabaseConnectionError(StorageError):
    def __init__(self, e):
        logger.error(e)


class DatabaseTableCreationError(StorageError):
    def __init__(self, e):
        logger.error(e)


class DatabaseQueryError(StorageError):
    def __init__(self, e):
        logger.error(e)


class DatabaseViewCreationError(StorageError):
    def __init__(self, e):
        logger.error(e)


class RTCoreError(Exception):
    def __init__(self, e):
        logger.error(e)


class OptionError(Exception):
    def __init__(self, e):
        logger.error(e)


class FileParsingError(Exception):
    def __init__(self, e):
        logger.error(e)


class WriteToStorageError(Exception):
    def __init__(self, e):
        logger.error(e)


class MultiprocessingError(Exception):
    def __init__(self, e):
        logger.error(e)


class ResultsProcessingError(Exception):
    def __init__(self, e):
        logger.error(e)


class OutputError(Exception):
    def __init__(self, e):
        logger.error(e)
