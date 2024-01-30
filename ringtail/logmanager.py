#   This module will handle any operational logging in Ringtail

import logging
import inspect
import traceback



class RTLogger:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(RTLogger, cls).__new__(cls)
            # Put any initialization here.
            cls.initialization(cls)
        return cls._instance
    
    def initialization(self, level = "WARNING", path = "rt_process_log.txt"):
        self.logger = logging.getLogger("ringtail")
        self.logger.setLevel(level)
        self.fileHandler = logging.FileHandler(path)
        self.logger.addHandler(self.fileHandler) 
        self.streamHandler = logging.StreamHandler()
        self.streamHandler.setLevel("WARNING")
        self.logger.addHandler(self.streamHandler) 
        self.streamFmt = logging.Formatter("%(levelname)-10s %(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        self.header(self, header="Starting a new Ringtail Process")
        self.fileHandler.setFormatter(self.streamFmt)

    def setLevel(self, level: str):
        self.logger.setLevel(level.upper())
        self.logger.debug("Log level changed to " + str(level))

    def level(self):
        levels ={10: "DEBUG",
                 20: "INFO",
                 30: "WARNING",
                 40: "ERROR",
                 50: "CRITICAL"}
        return levels[self.logger.level]

    def header(self, header):
        headFmt = logging.Formatter("%(asctime)30s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        self.fileHandler.setFormatter(headFmt)
        self.logger.critical("--------- " + header + " ---------")
        self.fileHandler.setFormatter(self.streamFmt)

    def debug(self, message):
        stacktrace = self.st_formatted(inspect.stack())
        self.logger.debug(stacktrace + ": " + message)

    def info(self, message):
        stacktrace = self.st_formatted(inspect.stack())
        self.logger.info(stacktrace + ": " + message)

    def warning(self, message):
        stacktrace = self.st_formatted(inspect.stack())
        self.logger.warning(stacktrace + ": " + message)

    def error(self, message):
        stacktrace = self.st_formatted(inspect.stack())
        self.logger.error(stacktrace + ": " + message)

    def critical(self, message):
        stacktrace = self.st_formatted(inspect.stack())
        self.logger.critical(stacktrace + ": " + message)

    def st_formatted(self, stack):
        stacktrace=""
        
    def st_formatted(self, stack):
        stacktrace=""
        for i in stack:
            fn = i.filename
            if fn.startswith("<frozen"):
                pass
            else:
                fn = fn.rsplit("/", 1)[1]
                lineno = i.lineno
                if fn == "exceptions.py" or fn =="logmanager.py":
                    pass
                else:
                    stacktrace += fn + "[" + str(lineno) + "]:"
        return stacktrace 


logger = RTLogger()
