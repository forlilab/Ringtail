#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail custom exceptions
#


class StorageError(Exception):
    pass


class DatabaseInsertionError(StorageError):
    pass


class DatabaseConnectionError(StorageError):
    pass


class DatabaseTableCreationError(StorageError):
    pass


class DatabaseQueryError(StorageError):
    pass


class DatabaseViewCreationError(StorageError):
    pass


class RTCoreError(Exception):
    pass


class OptionError(Exception):
    pass


class FileParsingError(Exception):
    pass


class WriteToStorageError(Exception):
    pass


class MultiprocessingError(Exception):
    pass


class ResultsProcessingError(Exception):
    pass


class OutputError(Exception):
    pass
