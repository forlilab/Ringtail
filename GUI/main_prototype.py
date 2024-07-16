import os
import sys
import util as u
from ringtail_prototype_ui import Ringtail_Prototype_UI
from PyQt6 import QtWidgets, QtGui
from PyQt6.QtWidgets import QApplication, QMainWindow
from ringtail import RingtailCore, RaccoonLogger
import resources_rc
from ui_log_handler import ConsoleWindowLogHandler
from file_browser import FileBrowser


class UI_MainWindow(Ringtail_Prototype_UI):
    def __init__(self):

        # initialize main window with the properties of the classes it inherits from (by calling super)
        super(UI_MainWindow, self).__init__()
        self.files = None
        self.filelists = None
        self.directories = None
        self.logger = RaccoonLogger(log_level_console="ERROR")

    def connectUI(self, MainWindow):
        ### set up the superficial UI (this method is inherited from ringtail_protptype_ui)
        self.setupUi(MainWindow)

        ### Connecting text boxes and buttons

        # region initWidget
        self.initRTButton.clicked.connect(self.pressedInitRT)
        self.selectPathsButton.setEnabled(False)
        self.submitResultsButton.setEnabled(False)
        consoleHandler = ConsoleWindowLogHandler()
        consoleHandler.sigLog.connect(self.logOutputTextBrowser.append)
        self.logger.logger.addHandler(consoleHandler)
        # endregion

        # region tabWidget

        # region resultsAdd
        self.selectPathsButton.clicked.connect(self.selectDockingResultPaths)
        self.submitResultsButton.clicked.connect(self.addDockingResults)
        # endregion

        # region databaseView

        # endregion

        # endregion

    def pressedInitRT(self):
        """
        Initializes the Ringtail core
        """
        if len(self.logFile.toPlainText()) > 0:
            log_file = self.logFile.toPlainText()
        else:
            log_file = None

        self.rtc = RingtailCore(
            db_file=self.dbFile.toPlainText(),
            storage_type=self.dbTypeDropdown.currentText(),
            docking_mode=self.dockingMode(),
            logging_level=self.logLevelDropdown.currentText(),
            logging_file=log_file,
        )
        # enable picking results paths once initialized
        self.selectPathsButton.setEnabled(True)
        # disable editing log file path after init
        self.logFile.setReadOnly(True)

    def dockingMode(self):
        """
        Determine docking mode from radio button

        Returns:
            str: docking mode
        """
        if self.adgpuButton.isChecked():
            self.filepattern = "*.dlg*"
            return "dlg"
        elif self.vinaButton.isChecked():
            self.filepattern = "*.pdbqt*"
            return "vina"

    def selectDockingResultPaths(self):
        # should open up a new window with file picking
        window = FileBrowser(self)
        window.show()

    def addDockingResults(self):
        self.rtc.add_results_from_files(
            file=u.QListWidget_to_list(self.files),
            file_path=u.QListWidget_to_list(self.directories),
            file_list=u.QListWidget_to_list(self.filelists),
        )
        self.files = None
        self.directories = None
        self.filelists = None


if __name__ == "__main__":
    import sys

    # Create a QT app widget
    app = QApplication(sys.argv)
    # set the icon that will appear on the dock
    app.setWindowIcon(QtGui.QIcon(":ringtail_head"))
    # create the QT widget main window that will hold and display the app
    MainWindow = QMainWindow()
    # create a user interface class that will get all the settings
    ui = UI_MainWindow()
    # connect the main window that can be shown, with the ui class that holds all the settings
    ui.connectUI(MainWindow)
    # show the main window
    MainWindow.show()
    # start/execute the application
    sys.exit(app.exec())
