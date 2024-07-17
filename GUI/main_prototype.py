import sys
import util as u
from ringtail_prototype_ui import Ringtail_Prototype_UI
from PyQt6 import QtWidgets, QtGui, QtCore
from PyQt6.QtWidgets import QApplication, QMainWindow
from ringtail import RingtailCore, RaccoonLogger
import resources_rc
import sqlite3
from ui_log_handler import ConsoleWindowLogHandler
from file_browser import FileBrowser
from worker import Worker


class UI_MainWindow(Ringtail_Prototype_UI):
    def __init__(self):

        # initialize main window with the properties of the classes it inherits from (by calling super)
        super(UI_MainWindow, self).__init__()
        self.files = None
        self.filelists = None
        self.directories = None
        self.logger = RaccoonLogger(log_level_console="DEBUG")

    def connectUI(self, MainWindow):
        ### set up the superficial UI (inherited from ringtail_protptype_ui)
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

        # TODO I'd like to add a count of number of files added
        # TODO connect files added to progress bar somehow
        # TODO add count for failed files
        # endregion

        # region databaseView
        self.threadpool = QtCore.QThreadPool()
        self.dbViewTab.setDisabled(True)
        self.tableWidget = QtWidgets.QTableWidget(self.dbViewTab)
        self.tableWidget.setObjectName("tableWidget")
        # figure out where to put table based on the other buttons
        tableNameButtonX, tableNameButtonY, _, buttonHeight = (
            self.tableNameButton.geometry().getRect()
        )
        tableX = tableNameButtonX
        tableY = tableNameButtonY + buttonHeight * 2
        # need to get tab widget dimensions for the full table size
        _, _, tabWidgetWidth, tabWidgetHeight = self.tabWidget.geometry().getRect()
        tableWidth = tabWidgetWidth - 2 * tableX
        tableHeight = tabWidgetHeight - tableY - 2 * tableX
        self.tableWidget.setGeometry(
            QtCore.QRect(tableX, tableY, tableWidth, tableHeight)
        )
        self.tableWidget.setColumnCount(10)
        self.tableWidget.setRowCount(10)

        # use lambda expresssions to how the query is sent depending on what button is pushed
        self.tableNameButton.clicked.connect(
            lambda: self.workerStart(
                f"SELECT * FROM {self.tableNameTextBox.toPlainText()}"
            )
        )
        # TODO this can probably be better used as a dropdown of available table names
        self.SQLQueryButton.clicked.connect(
            lambda: self.workerStart(self.SQLQueryTextBox.toPlainText())
        )
        # endregion

        # region filtering
        # only enable if initialized database has results
        self.filterTab.setDisabled(True)
        # endregion

        # endregion

    def workerStart(self, query):

        worker = Worker(lambda: self.loadDataBase(query))
        self.threadpool.start(worker)

    def loadDataBase(self, query):

        self.conn = sqlite3.connect(self.rtc.db_file)
        try:
            cursor = self.conn.execute(query)
        except Exception as e:
            print("    THIS IS THE ERROR: ", e)
            # TODO error popup that prints to log and can be cleared
        # TODO redo this stackoverflow method to get column iterated over
        row_len = []
        for i in cursor:
            row_len.append(len(i))
        self.col_num = max(row_len)
        self.tableWidget.setRowCount(0)
        self.tableWidget.setColumnCount(int(self.col_num))
        # set table headers
        cursor = self.conn.execute(query)
        columnNames = list(map(lambda x: x[0], cursor.description))
        self.tableWidget.insertRow(0)
        for col, colName in enumerate(columnNames):
            self.tableWidget.setItem(0, col, QtWidgets.QTableWidgetItem(str(colName)))
        # actually fetches the data and inserts it to the table widget
        for row, row_data in enumerate(cursor):
            self.tableWidget.insertRow(row + 1)
            for col, col_data in enumerate(row_data):
                self.tableWidget.setItem(
                    row + 1, col, QtWidgets.QTableWidgetItem(str(col_data))
                )

        self.conn.close()

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
        self.dbViewTab.setEnabled(True)
        if self.dbNumOfResults() > 0:
            self.filterTab.setEnabled(True)
        else:
            self.filterTab.setDisabled(True)

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

        if self.dbNumOfResults() > 0:
            self.filterTab.setEnabled(True)
        else:
            self.filterTab.setDisabled(True)

    def getValueRange(self, columnName: str = "docking_score") -> tuple:
        """
        Method to get max and min value from a numerical column in the database

        Args:
            columnName (str): name of column to be queried

        Returns:
            tuple: (min, max) values of the selected column
        """
        with self.rtc.storageman:
            max, min = self.rtc.storageman._run_query(
                f"SELECT MIN({columnName}), MAX({columnName}) FROM Results"
            ).fetchall()[0]

        return (min, max)

    def dbNumOfResults(self):
        """
        Count number of results in the database

        Returns:
            int: number of results in the results table
        """
        with self.rtc.storageman:
            return self.rtc.storageman._run_query(
                "SELECT COUNT(*) FROM Results"
            ).fetchone()[0]


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
