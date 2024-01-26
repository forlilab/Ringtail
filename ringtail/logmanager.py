#   This module will handle any operational logging in Ringtail

import logging
import inspect



class RTLogger:
    def __init__(self, level = "WARNING", path = "rt_process_log.txt"):
        self.logger = logging.getLogger("ringtail")
        self.logger.setLevel(level)
        self.fileHandler = logging.FileHandler(path)
        self.logger.addHandler(self.fileHandler) 

    def setLevel(self, level: str):
        self.logger.setLevel(level.upper())

    def level(self):
        levels ={10: "DEBUG",
                 20: "INFO",
                 30: "WARNING",
                 40: "ERROR",
                 50: "CRITICAL"}
        return levels[self.logger.level]

    def initialize(self):
        initFmt = logging.Formatter("%(asctime)30s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        self.fileHandler.setFormatter(initFmt)
        self.logger.critical("--------- NEW RINGTAIL PROCESS ---------")
        streamFmt = logging.Formatter("%(levelname)-10s %(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        self.fileHandler.setFormatter(streamFmt)

    def debug(self, message: str):
        stacktrace = str(inspect.stack()[1].filename) + ":" + str(inspect.stack()[1].function)+"["+str(inspect.stack()[1].lineno) + "]"
        self.logger.debug(stacktrace + " >>> " + message)

    def info(self, message):
        stacktrace = str(inspect.stack()[1].filename) + ":" + str(inspect.stack()[1].function)+"["+str(inspect.stack()[1].lineno) + "]"
        self.logger.info(stacktrace + " >>> " + message)

    def warning(self, message):
        stacktrace = str(inspect.stack()[1].filename) + ":" + str(inspect.stack()[1].function)+"["+str(inspect.stack()[1].lineno) + "]"
        self.logger.warning(stacktrace + " >>> " + message)

    def error(self, message):
        stacktrace = str(inspect.stack()[1].filename) + ":" + str(inspect.stack()[1].function)+"["+str(inspect.stack()[1].lineno) + "]"
        self.logger.error(stacktrace + " >>> " + message)

    def critical(self, message):
        stacktrace = str(inspect.stack()[1].filename) + ":" + str(inspect.stack()[1].function)+"["+str(inspect.stack()[1].lineno) + "]"
        self.logger.critical(stacktrace + " >>> " + message)

logger = RTLogger()
