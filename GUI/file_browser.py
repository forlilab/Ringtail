import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog
from PyQt5.uic import loadUi


class MainWindow(QDialog):
    def __init__(self):
        super(MainWindow, self).__init__()
        loadUi("GUI/file_dialog.ui", self)
        self.browse_files_button.clicked.connect(self.browsefiles)
        self.browse_directories_button.clicked.connect(self.handleChooseDirectories)
        self.selected_data = []

        # self.button = QtWidgets.QPushButton("Choose Directories")
        # self.button.clicked.connect(self.handleChooseDirectories)
        # self.listWidget = QtWidgets.QListWidget()
        # layout = QtWidgets.QVBoxLayout(self)
        # layout.addWidget(self.listWidget)
        # layout.addWidget(self.button)

    def browsefiles(self):

        filenames, _ = QFileDialog.getOpenFileNames(
            self,
            "QFileDialog.getOpenFileNames()",
            "",
        )
        if filenames:
            self.selected_data.extend(filenames)
            # add latest items on top of list
            # but this should always reflect self.selected_data I feel, and update selected data if the list view is changed, so probably a better approach needed in future
            self.data_list.insertItems(0, filenames)
            self.file_browser_text.setText(filenames[0])

    def browsedirectories(self):

        directory = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        print(directory)
        self.selected_data.extend(directory)
        self.data_list.insertItems(0, directory)
        self.file_browser_text.setText(directory)

    def handleChooseDirectories(self):
        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle("Choose Directories")
        dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
        dialog.setFileMode(QtWidgets.QFileDialog.DirectoryOnly)
        for view in dialog.findChildren((QtWidgets.QListView, QtWidgets.QTreeView)):
            if isinstance(view.model(), QtWidgets.QFileSystemModel):
                view.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)

        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            # self.listWidget.clear()
            # self.listWidget.addItems(dialog.selectedFiles())
            self.data_list.insertItems(0, dialog.selectedFiles())
            self.file_browser_text.setText(dialog.selectedFiles()[0])
        dialog.deleteLater()


app = QApplication(sys.argv)
mainwindow = MainWindow()
widget = QtWidgets.QStackedWidget()
widget.addWidget(mainwindow)
widget.setFixedWidth(400)
widget.setFixedHeight(300)
widget.show()
sys.exit(app.exec_())
