import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog
from PyQt5.uic import loadUi


class MainWindow(QDialog):
    def __init__(self):
        super(MainWindow, self).__init__()
        loadUi("GUI/file_dialog.ui", self)
        self.browse_files_button.clicked.connect(self.browsefiles)
        self.browse_directories_button.clicked.connect(self.browsedirectories)
        self.browse_filelists_button.clicked.connect(self.browsefilelists)
        self.delete_file_button.clicked.connect(self.delete_file_item)
        self.delete_directory_button.clicked.connect(self.delete_directory_item)
        self.delete_filelist_button.clicked.connect(self.delete_filelist_item)

    def selectfiles(self, extension: str = ""):
        filenames, _ = QFileDialog.getOpenFileNames(
            self,
            "QFileDialog.getOpenFileNames()",
            "",
            extension,
        )
        if filenames:
            # add latest items on top of list
            if "txt" in extension:
                self.filelist_list.insertItems(0, filenames)
            else:
                self.file_list.insertItems(0, filenames)
            self.file_browser_text.setText(filenames[0])

    def browsefiles(self):
        # this should be set by the chosen docking mode
        return self.selectfiles("*.dlg*")

    def browsedirectories(self):
        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle("Choose Directories")
        dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
        dialog.setFileMode(QtWidgets.QFileDialog.DirectoryOnly)
        for view in dialog.findChildren((QtWidgets.QListView, QtWidgets.QTreeView)):
            if isinstance(view.model(), QtWidgets.QFileSystemModel):
                view.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)

        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            self.directory_list.insertItems(0, dialog.selectedFiles())
            self.file_browser_text.setText(dialog.selectedFiles()[0])
        dialog.deleteLater()

    def browsefilelists(self):
        self.selectfiles("*.txt")

    def delete_file_item(self):
        # get index of the items, and remove from bottom up so index doesn't change
        indices = [x.row() for x in self.file_list.selectedIndexes()]
        # then remove each item from bottom up (so indices don't change)
        for index in reversed(indices):
            self.file_list.takeItem(index)

    def delete_directory_item(self):
        # get index of the items, and remove from bottom up so index doesn't change
        indices = [x.row() for x in self.directory_list.selectedIndexes()]

        # then remove each item from bottom up (so indices don't change)
        for index in reversed(indices):
            self.directory_list.takeItem(index)

    def delete_filelist_item(self):
        # get index of the items, and remove from bottom up so index doesn't change
        indices = [x.row() for x in self.filelist_list.selectedIndexes()]

        # then remove each item from bottom up (so indices don't change)
        for index in reversed(indices):
            self.filelist_list.takeItem(index)


app = QApplication(sys.argv)
mainwindow = MainWindow()
widget = QtWidgets.QStackedWidget()
widget.addWidget(mainwindow)
widget.setFixedWidth(800)
widget.setFixedHeight(600)
widget.show()
sys.exit(app.exec_())
