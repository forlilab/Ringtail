import sys
import logging
from PyQt5 import QtCore, QtWidgets

def setup_logger(debug, verbose):
    level = logging.DEBUG
    logging.basicConfig(
        level=level, stream=sys.stdout, filemode="w", format="%(message)s"
    )
    levels = {10: "debug", 20: "info", 30: "warning"}
    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.getLogger().setLevel(level)
    logging.info(f"Logging level set to {levels[level]}")


class QTextEditLogger(logging.Handler):
    def __init__(self, parent):
        super().__init__()
        self.widget = QtWidgets.QPlainTextEdit(parent)
        self.widget.setReadOnly(True)
        
    def emit(self, record):
        msg = self.format(record)
        self.widget.appendPlainText(msg)
        