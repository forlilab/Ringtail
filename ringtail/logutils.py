import logging
import os
from .util import caller_info


class RaccoonLogger:
    """This class implements the common functionalities for classes that need to
    provide logging facilities"""

    def __init__(self):
        """Initialize the logger"""
        self._owner = None
        self.setup_logger()

    def setup_logger(
        self,
        log_file: str = "RingtailLogger",
        log_console: bool = True,
        log_level: str = logging.WARNING,
        log_level_console: str = logging.WARNING,
        custom_logger_name: str = "RingtailLogger",
    ) -> None:
        """
        Setting up the logger.

        Args:
            log_file (str): name of log file, 'none' will not create a logfile
            log_console (bool): whether or not logger logs to console
            log_level (str): minimum logging level of the logger and file handler
            log_level_console (str): minimum logging level of the console, currently set with 'log_level'
            custom_logger_name (str): unique name of the logger
        """
        if self._owner is not None:
            return
        self._owner = caller_info()

        # access the logger; the logger module implements the named
        # singleton for logging, so if the logger is already available, it
        # will pass the existing one, or initialize it if necessary
        self.logger = logging.getLogger(custom_logger_name)
        # set the log level for the overal logger
        self.logger.setLevel(log_level)
        # configure the optional log file, if provided
        if self.logger.level == "DEBUG":
            self.add_filehandler(log_file)
        else:
            self._log_fp = None
        # initialize the console
        if log_console is True:
            self.log_console = logging.StreamHandler()
            self.log_console.setLevel(log_level_console)
            console_formatter = logging.Formatter("%(levelname)s - %(message)s")
            self.log_console.setFormatter(console_formatter)
            self.logger.addHandler(self.log_console)
        else:
            self.log_console = None

        self.logger.info(
            f'[ {custom_logger_name} ] log system succcessfully initialized. [ file log: "{log_file}" | console: {log_console}]'
        )

    def level(self):
        """
        Returns:
            str: current level of the logger
        """
        levels = {10: "DEBUG", 20: "INFO", 30: "WARNING", 40: "ERROR", 50: "CRITICAL"}
        return levels[self.logger.level]

    def add_filehandler(self, log_file: str = "ringtail", level: str = logging.DEBUG):
        """
        Will add file handler to an existing logging object and produce a '.log' file.

        Args:
            log_file (str): Name to include in the log file name, will prepend date and time.
            level (str): Log level of log file, defaults to DEBUG.
        """
        import datetime

        dt = datetime.datetime.now()
        filename = dt.strftime("%Y%m%d-%H%M%S") + "_" + log_file + ".log"
        self._log_fp = logging.FileHandler(filename)
        if level is not None:
            self._log_fp.setLevel(level)
        else:
            self._log_fp.setLevel(self.level())
        log_file_formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s"
        )
        self._log_fp.setFormatter(log_file_formatter)
        self.logger.addHandler(self._log_fp)

    def set_level(self, log_level: str):
        """Sets level of the logger, and if debug will print level change to log.

        Args:
            level (str): lowest level the logger will record, either DEBUG, INFO, or WARNING
        """
        if type(log_level) == str:
            log_level = log_level.upper()
        else:
            self.logger.warning(
                f"{log_level} is not a a string. Please use one of the following options: DEBUG, INFO, or WARNING. Logger level reverted to {self.level()}."
            )
            return
        if log_level not in ["DEBUG", "INFO", "WARNING"]:
            self.logger.warning(
                f"{log_level} is not a valid logging level option. Logger level reverted to {self.level()}."
            )
            return
        elif log_level != self.level():
            self.logger.setLevel(log_level)
            if self._log_fp is not None:
                self._log_fp.setLevel(log_level)
            if self.log_console is not None:
                self.log_console.setLevel(log_level)
            self.logger.debug("Log level changed to " + str(log_level))

    def get_caller(self, stack_level: int = 2):
        import inspect

        """
        Method to get basic information about the module and line no calling a certain function.

        Args:
            stack_level (int): what level in the stack you want. 0 is this method, 1 is its caller (typically the logutil.log method), 2 the caller of the caller (usually the package module)

        Returns:
            str: name of the module/file calling
            int: line number in the module/file calling
        """
        module_path = inspect.stack()[stack_level].filename
        module = os.path.basename(module_path)
        line = inspect.stack()[stack_level].lineno
        return module, line

    def debug(self, message):
        """
        Lowest level of log message

        Args:
            message (str): message to be logged
        """
        module, line = self.get_caller()
        self.logger.debug(module + ":" + str(line) + " - " + message)

    def info(self, message):
        """
        Second lowest level of log message

        Args:
            message (str): message to be logged
        """
        module, line = self.get_caller()
        self.logger.info(module + ":" + str(line) + " - " + message)

    def warning(self, message):
        """
        Medium level of log message

        Args:
            message (str): message to be logged
        """
        module, line = self.get_caller()
        self.logger.warning(module + ":" + str(line) + " - " + message)

    def error(self, message):
        """
        For logging errors

        Args:
            message (str): message to be logged. A stack trace will be appended automatically.
        """
        module, line = self.get_caller()
        message = module + ":" + str(line) + " - " + message
        self.logger.error(message + "\n", exc_info=1)

    def critical(self, message):
        """
        For logging critical errors and failures. A stack trace will be appended automatically.

        Args:
            message (str): message to be logged
        """
        module, line = self.get_caller()
        message = module + ":" + str(line) + " - " + message
        self.logger.critical(message + "\n", exc_info=1)


LOGGER = RaccoonLogger()
