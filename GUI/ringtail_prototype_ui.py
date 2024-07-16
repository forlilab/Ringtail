# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ringtail_prototype.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QMainWindow


class Ringtail_Prototype_UI(QMainWindow):
    def __init__(self):
        super(Ringtail_Prototype_UI, self).__init__()

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(881, 965)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.initwidget = QtWidgets.QWidget(self.centralwidget)
        self.initwidget.setGeometry(QtCore.QRect(10, 10, 861, 251))
        self.initwidget.setObjectName("initwidget")
        self.logOutputTextBrowser = QtWidgets.QTextBrowser(self.initwidget)
        self.logOutputTextBrowser.setGeometry(QtCore.QRect(340, 30, 491, 192))
        self.logOutputTextBrowser.setObjectName("logOutputTextBrowser")
        self.logOutputLabel = QtWidgets.QLabel(self.initwidget)
        self.logOutputLabel.setGeometry(QtCore.QRect(340, 0, 161, 16))
        self.logOutputLabel.setObjectName("logOutputLabel")
        self.logLevelLabel = QtWidgets.QLabel(self.initwidget)
        self.logLevelLabel.setGeometry(QtCore.QRect(190, 10, 91, 16))
        self.logLevelLabel.setObjectName("logLevelLabel")
        self.dbFile = QtWidgets.QPlainTextEdit(self.initwidget)
        self.dbFile.setGeometry(QtCore.QRect(10, 150, 191, 71))
        self.dbFile.setObjectName("dbFile")
        self.initRTButton = QtWidgets.QCommandLinkButton(self.initwidget)
        self.initRTButton.setGeometry(QtCore.QRect(197, 160, 121, 50))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(15)
        self.initRTButton.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(
            QtGui.QPixmap("../Ringtail_GUI_old/resources/ringtail_head.png"),
            QtGui.QIcon.Mode.Normal,
            QtGui.QIcon.State.Off,
        )
        self.initRTButton.setIcon(icon)
        self.initRTButton.setObjectName("initRTButton")
        self.dbTypeDropdown = QtWidgets.QComboBox(self.initwidget)
        self.dbTypeDropdown.setGeometry(QtCore.QRect(0, 30, 81, 22))
        self.dbTypeDropdown.setObjectName("dbTypeDropdown")
        self.dbTypeDropdown.addItem("")
        self.dbTypeDropdown.addItem("")
        self.dbTypeDropdown.addItem("")
        self.logFileLabel = QtWidgets.QLabel(self.initwidget)
        self.logFileLabel.setGeometry(QtCore.QRect(150, 70, 161, 20))
        self.logFileLabel.setObjectName("logFileLabel")
        self.vinaButton = QtWidgets.QRadioButton(self.initwidget)
        self.vinaButton.setGeometry(QtCore.QRect(0, 90, 91, 21))
        self.vinaButton.setObjectName("vinaButton")
        self.logFile = QtWidgets.QPlainTextEdit(self.initwidget)
        self.logFile.setGeometry(QtCore.QRect(150, 90, 161, 31))
        self.logFile.setPlainText("")
        self.logFile.setObjectName("logFile")
        self.dbTypeLabel = QtWidgets.QLabel(self.initwidget)
        self.dbTypeLabel.setGeometry(QtCore.QRect(0, 10, 91, 16))
        self.dbTypeLabel.setObjectName("dbTypeLabel")
        self.adgpuButton = QtWidgets.QRadioButton(self.initwidget)
        self.adgpuButton.setEnabled(True)
        self.adgpuButton.setGeometry(QtCore.QRect(0, 60, 151, 21))
        self.adgpuButton.setChecked(True)
        self.adgpuButton.setObjectName("adgpuButton")
        self.logLevelDropdown = QtWidgets.QComboBox(self.initwidget)
        self.logLevelDropdown.setGeometry(QtCore.QRect(190, 30, 81, 22))
        self.logLevelDropdown.setObjectName("logLevelDropdown")
        self.logLevelDropdown.addItem("")
        self.logLevelDropdown.addItem("")
        self.logLevelDropdown.addItem("")
        self.dbFileLabel = QtWidgets.QLabel(self.initwidget)
        self.dbFileLabel.setGeometry(QtCore.QRect(10, 130, 121, 16))
        self.dbFileLabel.setObjectName("dbFileLabel")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(10, 270, 851, 651))
        self.tabWidget.setObjectName("tabWidget")
        self.inputTab = QtWidgets.QWidget()
        self.inputTab.setObjectName("inputTab")
        self.numOfFiles = QtWidgets.QPlainTextEdit(self.inputTab)
        self.numOfFiles.setGeometry(QtCore.QRect(480, 40, 121, 41))
        self.numOfFiles.setObjectName("numOfFiles")
        self.numOfDirs = QtWidgets.QPlainTextEdit(self.inputTab)
        self.numOfDirs.setGeometry(QtCore.QRect(480, 100, 121, 41))
        self.numOfDirs.setObjectName("numOfDirs")
        self.numOfFilelists = QtWidgets.QPlainTextEdit(self.inputTab)
        self.numOfFilelists.setGeometry(QtCore.QRect(480, 160, 121, 41))
        self.numOfFilelists.setObjectName("numOfFilelists")
        self.numOfFilesLabel = QtWidgets.QLabel(self.inputTab)
        self.numOfFilesLabel.setGeometry(QtCore.QRect(620, 50, 181, 16))
        self.numOfFilesLabel.setObjectName("numOfFilesLabel")
        self.numOfDirsLabel = QtWidgets.QLabel(self.inputTab)
        self.numOfDirsLabel.setGeometry(QtCore.QRect(620, 110, 181, 16))
        self.numOfDirsLabel.setObjectName("numOfDirsLabel")
        self.numOfFilelistsLabel = QtWidgets.QLabel(self.inputTab)
        self.numOfFilelistsLabel.setGeometry(QtCore.QRect(620, 170, 181, 16))
        self.numOfFilelistsLabel.setObjectName("numOfFilelistsLabel")
        self.selectPathsButton = QtWidgets.QPushButton(self.inputTab)
        self.selectPathsButton.setGeometry(QtCore.QRect(100, 80, 251, 91))
        self.selectPathsButton.setObjectName("selectPathsButton")
        self.submitResultsButton = QtWidgets.QPushButton(self.inputTab)
        self.submitResultsButton.setGeometry(QtCore.QRect(100, 330, 251, 91))
        self.submitResultsButton.setObjectName("submitResultsButton")
        self.resultsAddProgressBar = QtWidgets.QProgressBar(self.inputTab)
        self.resultsAddProgressBar.setGeometry(QtCore.QRect(130, 480, 191, 23))
        self.resultsAddProgressBar.setProperty("value", 24)
        self.resultsAddProgressBar.setObjectName("resultsAddProgressBar")
        self.numOfFailedLabel = QtWidgets.QLabel(self.inputTab)
        self.numOfFailedLabel.setGeometry(QtCore.QRect(620, 420, 181, 16))
        self.numOfFailedLabel.setObjectName("numOfFailedLabel")
        self.numOfFailed = QtWidgets.QPlainTextEdit(self.inputTab)
        self.numOfFailed.setGeometry(QtCore.QRect(480, 410, 121, 41))
        self.numOfFailed.setObjectName("numOfFailed")
        self.numOfResultsLabel = QtWidgets.QLabel(self.inputTab)
        self.numOfResultsLabel.setGeometry(QtCore.QRect(620, 360, 181, 16))
        self.numOfResultsLabel.setObjectName("numOfResultsLabel")
        self.numOfResults = QtWidgets.QPlainTextEdit(self.inputTab)
        self.numOfResults.setGeometry(QtCore.QRect(480, 350, 121, 41))
        self.numOfResults.setObjectName("numOfResults")
        self.tabWidget.addTab(self.inputTab, "")
        self.dbViewTab = QtWidgets.QWidget()
        self.dbViewTab.setObjectName("dbViewTab")
        self.tableNameButton = QtWidgets.QPushButton(self.dbViewTab)
        self.tableNameButton.setGeometry(QtCore.QRect(10, 110, 120, 26))
        self.tableNameButton.setObjectName("tableNameButton")
        self.tableNameTextBox = QtWidgets.QPlainTextEdit(self.dbViewTab)
        self.tableNameTextBox.setGeometry(QtCore.QRect(10, 60, 121, 41))
        self.tableNameTextBox.setObjectName("tableNameTextBox")
        self.tableNameLabel = QtWidgets.QLabel(self.dbViewTab)
        self.tableNameLabel.setGeometry(QtCore.QRect(10, 30, 91, 16))
        self.tableNameLabel.setObjectName("tableNameLabel")
        self.SQLQueryLabel = QtWidgets.QLabel(self.dbViewTab)
        self.SQLQueryLabel.setGeometry(QtCore.QRect(182, 30, 151, 16))
        self.SQLQueryLabel.setObjectName("SQLQueryLabel")
        self.SQLQueryTextBox = QtWidgets.QPlainTextEdit(self.dbViewTab)
        self.SQLQueryTextBox.setGeometry(QtCore.QRect(179, 60, 661, 41))
        self.SQLQueryTextBox.setObjectName("SQLQueryTextBox")
        self.orLabel = QtWidgets.QLabel(self.dbViewTab)
        self.orLabel.setGeometry(QtCore.QRect(142, 30, 30, 16))
        self.orLabel.setObjectName("orLabel")
        self.SQLQueryButton = QtWidgets.QPushButton(self.dbViewTab)
        self.SQLQueryButton.setGeometry(QtCore.QRect(180, 110, 658, 26))
        self.SQLQueryButton.setObjectName("SQLQueryButton")
        self.dbViewWidget = QtWidgets.QWidget(self.dbViewTab)
        self.dbViewWidget.setGeometry(QtCore.QRect(10, 150, 821, 461))
        self.dbViewWidget.setObjectName("dbViewWidget")
        self.dbViewLabel = QtWidgets.QLabel(self.dbViewWidget)
        self.dbViewLabel.setGeometry(QtCore.QRect(350, 180, 58, 16))
        self.dbViewLabel.setObjectName("dbViewLabel")
        self.tabWidget.addTab(self.dbViewTab, "")
        self.plotTab = QtWidgets.QWidget()
        self.plotTab.setObjectName("plotTab")
        self.graphicsView = QtWidgets.QGraphicsView(self.plotTab)
        self.graphicsView.setGeometry(QtCore.QRect(430, 41, 371, 351))
        self.graphicsView.setObjectName("graphicsView")
        self.updatePlotButton = QtWidgets.QCommandLinkButton(self.plotTab)
        self.updatePlotButton.setGeometry(QtCore.QRect(10, 110, 168, 41))
        self.updatePlotButton.setObjectName("updatePlotButton")
        self.plotAllDataButton = QtWidgets.QRadioButton(self.plotTab)
        self.plotAllDataButton.setEnabled(True)
        self.plotAllDataButton.setGeometry(QtCore.QRect(10, 30, 151, 21))
        self.plotAllDataButton.setChecked(True)
        self.plotAllDataButton.setObjectName("plotAllDataButton")
        self.plot_data_button_group = QtWidgets.QButtonGroup(MainWindow)
        self.plot_data_button_group.setObjectName("plot_data_button_group")
        self.plot_data_button_group.addButton(self.plotAllDataButton)
        self.plotFilterDataButton = QtWidgets.QRadioButton(self.plotTab)
        self.plotFilterDataButton.setGeometry(QtCore.QRect(10, 60, 91, 21))
        self.plotFilterDataButton.setObjectName("plotFilterDataButton")
        self.plot_data_button_group.addButton(self.plotFilterDataButton)
        self.tabWidget.addTab(self.plotTab, "")
        self.filterTab = QtWidgets.QWidget()
        self.filterTab.setObjectName("filterTab")
        self.tabWidget.addTab(self.filterTab, "")
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.logOutputLabel.setText(_translate("MainWindow", "Log messages"))
        self.logLevelLabel.setText(_translate("MainWindow", "Logging level"))
        self.dbFile.setPlainText(_translate("MainWindow", "output.db"))
        self.initRTButton.setText(
            _translate("MainWindow", "Start your \n" "Ringtail Core")
        )
        self.dbTypeDropdown.setItemText(0, _translate("MainWindow", "sqlite"))
        self.dbTypeDropdown.setItemText(1, _translate("MainWindow", "sqlalchemy"))
        self.dbTypeDropdown.setItemText(2, _translate("MainWindow", "postgres"))
        self.logFileLabel.setText(
            _translate("MainWindow", "Logging file name (optional)")
        )
        self.vinaButton.setText(_translate("MainWindow", "vina (pdbqt)"))
        self.dbTypeLabel.setText(_translate("MainWindow", "Database type"))
        self.adgpuButton.setText(_translate("MainWindow", "AutoDock GPU (dlg)"))
        self.logLevelDropdown.setItemText(0, _translate("MainWindow", "warning"))
        self.logLevelDropdown.setItemText(1, _translate("MainWindow", "info"))
        self.logLevelDropdown.setItemText(2, _translate("MainWindow", "debug"))
        self.dbFileLabel.setText(_translate("MainWindow", "Database path"))
        self.numOfFilesLabel.setText(_translate("MainWindow", "files chosen"))
        self.numOfDirsLabel.setText(_translate("MainWindow", "directories chosen"))
        self.numOfFilelistsLabel.setText(_translate("MainWindow", "file lists chosen"))
        self.selectPathsButton.setText(_translate("MainWindow", "Select results"))
        self.submitResultsButton.setText(
            _translate("MainWindow", "Submit results to database")
        )
        self.numOfFailedLabel.setText(_translate("MainWindow", "failed files"))
        self.numOfResultsLabel.setText(_translate("MainWindow", "results added"))
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.inputTab),
            _translate("MainWindow", "Add docking results"),
        )
        self.tableNameButton.setText(_translate("MainWindow", "Show table"))
        self.tableNameLabel.setText(_translate("MainWindow", "Table name"))
        self.SQLQueryLabel.setText(_translate("MainWindow", "SQL formatted query"))
        self.orLabel.setText(_translate("MainWindow", "OR"))
        self.SQLQueryButton.setText(_translate("MainWindow", "Show results of query"))
        self.dbViewLabel.setText(_translate("MainWindow", "(db view)"))
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.dbViewTab),
            _translate("MainWindow", "View database"),
        )
        self.updatePlotButton.setText(_translate("MainWindow", "Update plot"))
        self.plotAllDataButton.setText(_translate("MainWindow", "All data"))
        self.plotFilterDataButton.setText(_translate("MainWindow", "Filtered data"))
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.plotTab), _translate("MainWindow", "Plotting")
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.filterTab),
            _translate("MainWindow", "Filtering"),
        )
