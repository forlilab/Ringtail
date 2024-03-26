import logging
import logging.handlers
import inspect
import datetime

class RTLogger:
    """
    RTLogger is a global singleton class for code in Ringtail. A new header is written 
    each time a new ringtail sore object is started. It is based on the python logging 
    library. Different log levels have different treatment, e.g., error level logs includes
    a complete and readable stack trace. 

    The log is written to "rt_process_log.txt", and the file will be saved to "Ringtail/logfiles/".

    """

    _instance = None
    def __new__(cls):
        #Ensures class is only instantiated once
        if cls._instance is None:
            cls._instance = super(RTLogger, cls).__new__(cls)
            # Put any initialization here.
            dt = datetime.datetime.now()
            filename = dt.strftime("%Y%m%d-%H%M%S") + "_ringtail-process-log.txt"
            cls.initialization(cls, filename=filename)
        return cls._instance
    
    #NOTE default log level set to debug during testing
    def initialization(self, level = "DEBUG", filename = "rt_process_log.txt"):
        #TODO get log level default from options
        """ 
        Options for instantiation of the logger. 
        """
        self.logger = logging.getLogger("ringtail")
        self.logger.setLevel(level)
        self.fileHandler = logging.FileHandler(filename)
        self.logger.addHandler(self.fileHandler) 
        self.streamHandler = logging.StreamHandler()
        self.streamHandler.setLevel("WARNING")
        self.logger.addHandler(self.streamHandler) 
        self.streamFmt = logging.Formatter("%(levelname)-10s %(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        self.header(self, header="Starting a new Ringtail Process")
        self.fileHandler.setFormatter(self.streamFmt)

    def setLevel(self, level: str):
        """ Sets level of the logger and prints to log (if debug)."""
        if type(level) == str:
            level = level.upper()
        else:
            logger.warning(f"{level} is not a a string. Please use one of the following options: debug, info, or warning. Logger level reverted to {self.level()}.")
            return
        
        if level.upper() not in ["DEBUG", "INFO", "WARNING","ERROR","CRITICAL"]:
            logger.warning(f"{level} is not a valid logging level option. Logger level reverted to {self.level()}.")
            return
        elif level != self.logger.level:
            self.logger.setLevel(level.upper())
            self.logger.debug("Log level changed to " + str(level))

    def level(self):
        levels ={10: "DEBUG",
                 20: "INFO",
                 30: "WARNING",
                 40: "ERROR",
                 50: "CRITICAL"}
        return levels[self.logger.level]

    def filename(self):
        return self.fileHandler.baseFilename
    
    def header(self, header):
        headFmt = logging.Formatter("%(asctime)30s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        self.fileHandler.setFormatter(headFmt)
        self.logger.critical("--------- " + header + " ---------")
        self.fileHandler.setFormatter(self.streamFmt)

    def debug(self, message):
        self.logger.debug(message)

    def info(self, message):
        self.logger.info(message)

    def warning(self, message):
        self.logger.warning(message)

    def error(self, message):
        stacktrace = self.st_formatted(inspect.stack())
        self.logger.error(stacktrace + ": " + message)

    def critical(self, message):
        stacktrace = self.st_formatted(inspect.stack())
        self.logger.critical(stacktrace + ": " + message)
        
    def st_formatted(self, stack):
        """
        Method to format python FrameInfo object to this format:
        file[lineno]:file[lineno]file[lineno]: (first file is last in the stack)
        """
        stacktrace=""
        for i in stack:
            
            fn: str = i.filename
            if fn.startswith("<frozen"):
                pass
            elif "site-packages" in fn:
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
