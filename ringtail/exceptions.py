#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail custom exceptions
#

class DatabaseError(Exception):
    pass


class DatabaseInsertionError(DatabaseError):
    pass


class DatabaseConnectionError(DatabaseError):
    pass


class DatabaseTableCreationError(DatabaseError):
    pass


class DatabaseQueryError(DatabaseError):
    pass


class DatabaseViewCreationError(DatabaseError):
    pass


class VirtualScreeningError(Exception):
    pass


class OptionError(Exception):
    pass


class FileParsingError(Exception):
    pass


class WriteToDatabaseError(Exception):
    pass


class MultiprocessingError(Exception):
    pass


class ResultsProcessingError(Exception):
    pass


class OutputError(Exception):
    pass
