#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail custom exceptions
#


class RingtailError(Exception):
    """Base class for ringtail exceptions"""


class StorageError(RingtailError):
    pass


class MergeError(StorageError):
    pass


class VersionError(StorageError):
    pass


class DatabaseInsertionError(StorageError):
    pass


class DatabaseConnectionError(StorageError):
    pass


class DatabaseQueryError(StorageError):
    pass


class RTCoreError(RingtailError):
    pass


class OptionError(RingtailError):
    pass


class FileParsingError(RingtailError):
    pass


class FileParsingErrorAdgpu(FileParsingError):
    pass


class FileParsingErrorPdbqt(FileParsingError):
    pass


class FileParsingErrorSdf(FileParsingError):
    pass


class WriteToStorageError(RingtailError):
    pass


class MultiprocessingError(RingtailError):
    pass


class ResultsProcessingError(RingtailError):
    pass


class OutputError(RingtailError):
    pass


class InteractionError(RingtailError):
    pass
