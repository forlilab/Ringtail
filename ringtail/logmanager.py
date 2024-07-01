#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail logging manager
#

import logging
import logging.handlers
import inspect
import datetime

from typing import Union

# https://gist.github.com/lee-pai-long/d3004225e1847b84acb4fbba0c2aea91
"""
            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
                    Version 2, December 2004

 Copyright (C) 2004 Sam Hocevar <sam@hocevar.net>

 Everyone is permitted to copy and distribute verbatim or modified
 copies of this license document, and changing it is allowed as long
 as the name is changed.

            DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. You just DO WHAT THE FUCK YOU WANT TO.
"""
import inspect


def caller_info(skip=1):
    """Get the name of a caller in the format module.class.method.

    https://gist.github.com/lee-pai-long/d3004225e1847b84acb4fbba0c2aea91
    Copied from: https://gist.github.com/techtonik/2151727

    :arguments:
        - skip (integer): Specifies how many levels of stack
                          to skip while getting caller name.
                          skip=1 means "who calls me",
                          skip=2 "who calls my caller" etc.

    :returns:
        - package (string): caller package.
        - module (string): caller module.
        - klass (string): caller classname if one otherwise None.
        - caller (string): caller function or method (if a class exist).
        - line (int): the line of the call.
        - An empty string is returned if skipped levels exceed stack height.
    """
    stack = inspect.stack()
    start = 0 + skip
    if len(stack) < start + 1:
        return ""
    parentframe = stack[start][0]

    # module and packagename.
    module_info = inspect.getmodule(parentframe)
    if module_info:
        mod = module_info.__name__.split(".")
        package = mod[0]
        try:
            module = mod[1]
        except:
            module = ""

    # class name.
    klass = None
    if "self" in parentframe.f_locals:
        klass = parentframe.f_locals["self"].__class__.__name__

    # method or function name.
    caller = None
    if parentframe.f_code.co_name != "<module>":  # top level usually
        caller = parentframe.f_code.co_name

    # call line.
    line = parentframe.f_lineno

    # Remove reference to frame
    # See: https://docs.python.org/3/library/inspect.html#the-interpreter-stack
    del parentframe

    return package, module, klass, caller, line


class RingtailLogger:
    """
    RTLogger is a global singleton class for code in Ringtail. A new header is written
    each time a new ringtail sore object is started. It is based on the python logging
    library. Different log levels have different treatment, e.g., error level logs includes
    a complete and readable stack trace.

    The log is written to "YYYYMMDD-hhmmss_ringtail-process-log.txt", and the file will be saved to current working directory.
    """

    _instance = {}

    def __init__(self):
        """initializing the logger"""
        self._owner = None
        self.logger_setup(log_file="ringtail-process-log")

    def logger_setup(
        self,
        log_file: Union[str, None] = None,
        log_console: bool = False,
        log_level=logging.DEBUG,
        log_level_console=logging.WARNING,
        custom_logger_name: str = "RingtailLogger",
    ) -> None:
        if not self._owner is None:
            print("[already initialized, skipping]")
            print("OWNER:>", self._owner)
            return
        self._owner = caller_info()
        print(
            "SETTING UP!",
            log_file,
            log_console,
            log_level,
            log_level_console,
            custom_logger_name,
        )
        print("CALLER>", self._owner)

        # access the logger; the logger module implements the named
        # singleton for logging, so if the logger is already available, it
        # will pass the existing one, or initialize it if necessary
        self.logger = logging.getLogger(custom_logger_name)
        # set the log level for the overal logger
        self.logger.setLevel(log_level)

        # configure the optional log file, if provided
        if log_file is not None:
            self._log_fp = logging.FileHandler(log_file)
            self._log_fp.setLevel(log_level)
            log_file_formatter = logging.Formatter(
                "%(asctime)s - %(levelname)s - %(message)s"
            )
            self._log_fp.setFormatter(log_file_formatter)
            self.logger.addHandler(self._log_fp)
        else:
            self._log_fp = None
        # initialize the console
        if log_console is True:
            self.log_console = logging.StreamHandler()
            self.log_console.setLevel(log_level_console)
            console_formatter = logging.Formatter(
                "%(levelname)s - %(filename)s - Line: %(lineno)d - %(message)s"
            )
            self.log_console.setFormatter(console_formatter)
            self.logger.addHandler(self.log_console)
        else:
            self.log_console = None
        self.logger.info(
            f'[ RingtailLogger ] log system succcessfully initialized. [ file log: "{log_file}" | console {log_console}]'
        )

    def log(self, message):
        print("LOG MESSAGE USING LOGGER ID", id(self))
        # self.logger.info(message)


logger = RingtailLogger()

#################

#     self, level: str = "WARNING", filename="ringtail-process-log"):
#     """Initializing and starting the logger.
#     Starts a file handler and a stream handler that prints to stdout

#     Args:
#         level (str): logger level
#         filename (str): filename of the log file
#     """
#     if self._active is True:
#         print("Logger already initialized, not performing additional initializations.")
#         return
#     print("CALLED ", filename)
#     self.logger = logging.getLogger()
#         self._instance = super(RTLogger, self).__new__(self)
#         dt = datetime.datetime.now()
#         filename = dt.strftime("%Y%m%d-%H%M%S") + "_" + filename + ".txt"
#         self.initialization(self, filename=filename)

#     self.logger = logging.getLogger("ringtail")
#     self.logger.setLevel(level)
#     self.fileHandler = logging.FileHandler(filename)
#     self.logger.addHandler(self.fileHandler)
#     self.streamHandler = logging.StreamHandler()
#     self.streamHandler.setLevel(level)
#     self.logger.addHandler(self.streamHandler)
#     self.streamFmt = logging.Formatter("%(levelname)-10s %(message)s")
#     self.fileFmt = logging.Formatter(
#         "%(levelname)-10s %(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
#     )
#     self.header(self, header="Starting a new Ringtail Process")
#     self.fileHandler.setFormatter(self.fileFmt)
#     self.streamHandler.setFormatter(self.streamFmt)

#     return self._instance

# def __new__(cls):
#     """Method to ensure singleton logger"""
#     if cls._instance is None:
#         cls._instance = super(RTLogger, cls).__new__(cls)
#         dt = datetime.datetime.now()
# filename = dt.strftime("%Y%m%d-%H%M%S") + "_ringtail-process-log.txt"
#         cls.initialization(cls, filename=filename)
#     return cls._instance

# def initialization(self, level="WARNING", filename="ringtail-process-log.txt"):
#     """
#     Options for instantiation of the logger.
#     Starts a file handler and a stream handler that prints to stdout

#     Args:
#         level (str): logger level
#         filename (str): filename of the log file
#     """
#     self.logger = logging.getLogger("ringtail")
#     self.logger.setLevel(level)
#     self.fileHandler = logging.FileHandler(filename)
#     self.logger.addHandler(self.fileHandler)
#     self.streamHandler = logging.StreamHandler()
#     self.streamHandler.setLevel(level)
#     self.logger.addHandler(self.streamHandler)
#     self.streamFmt = logging.Formatter("%(levelname)-10s %(message)s")
#     self.fileFmt = logging.Formatter(
#         "%(levelname)-10s %(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
#     )
#     self.header(self, header="Starting a new Ringtail Process")
#     self.fileHandler.setFormatter(self.fileFmt)
#     self.streamHandler.setFormatter(self.streamFmt)


# def header(self, header):
#     """
#     Formats header for the log file

#     Args:
#         header (str): header string to get formatted
#     """
#     headFmt = logging.Formatter(
#         "%(asctime)30s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
#     )
#     self.fileHandler.setFormatter(headFmt)
#     self.logger.critical("--------- " + header + " ---------")
#     self.fileHandler.setFormatter(self.fileFmt)


# def st_formatted(self, stack):
#     """
#     Method to format python FrameInfo object to this format:
#     file[lineno]:file[lineno]file[lineno]: (first file is last in the stack)

#     Args:
#         stack (traceback): stack to be formatted to custom traceback
#     """
#     stacktrace = ""
#     for i in stack:

#         fn: str = i.filename
#         if fn.startswith("<frozen"):
#             pass
#         elif "site-packages" in fn:
#             pass
#         else:
#             fn = fn.rsplit("/", 1)[1]
#             lineno = i.lineno
#             if fn == "exceptions.py" or fn == "logmanager.py":
#                 pass
#             else:
#                 stacktrace += fn + "[" + str(lineno) + "]:"
#     return stacktrace
