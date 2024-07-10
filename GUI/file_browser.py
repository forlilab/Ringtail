import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog
from PyQt5.uic import loadUi


class FileBrowser(QDialog):
    def __init__(self):
        super(FileBrowser, self).__init__()
        loadUi("GUI/file_dialog.ui", self)
        self.browse_files_button.clicked.connect(self.browsefiles)
        self.browse_directories_button.clicked.connect(self.browsedirectories)
        self.browse_filelists_button.clicked.connect(self.browsefilelists)
        self.delete_file_button.clicked.connect(
            lambda: self.delete_item(self.file_list)
        )
        self.delete_directory_button.clicked.connect(
            lambda: self.delete_item(self.directory_list)
        )
        self.delete_filelist_button.clicked.connect(
            lambda: self.delete_item(self.filelist_list)
        )

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

    def delete_item(self, list_widget):
        # get index of the items, and remove from bottom up so index doesn't change
        indices = [x.row() for x in list_widget.selectedIndexes()]
        # then remove each item from bottom up (so indices don't change)
        for index in reversed(indices):
            list_widget.takeItem(index)


if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    filebrowser = FileBrowser()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(filebrowser)
    widget.setFixedWidth(800)
    widget.setFixedHeight(600)
    widget.show()
    sys.exit(app.exec_())
