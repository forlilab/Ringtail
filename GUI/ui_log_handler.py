import logging
from PyQt6.QtCore import QObject, pyqtSignal


class ConsoleWindowLogHandler(logging.Handler, QObject):
    """
    Class to handle a UI logging handler that emits a signal to a text box
    """

    sigLog = pyqtSignal(str)

    def __init__(self):
        logging.Handler.__init__(self)
        QObject.__init__(self)

    def emit(self, logRecord):
        message = str(logRecord.getMessage())
        self.sigLog.emit(message)
