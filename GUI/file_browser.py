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
        self.selected_files = []
        self.selected_directories = []
        self.selected_filelists = []

    def selectfiles(self, extension: str = ""):
        filenames, _ = QFileDialog.getOpenFileNames(
            self,
            "QFileDialog.getOpenFileNames()",
            "",
            extension,
        )
        if filenames:

            # add latest items on top of list
            # but this should always reflect self.selected_data I feel, and update selected data if the list view is changed, so probably a better approach needed in future
            if "txt" in extension:
                self.filelist_list.insertItems(0, filenames)
                self.selected_filelists.extend(filenames)
            else:
                self.file_list.insertItems(0, filenames)
                self.selected_files.extend(filenames)
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
            self.selected_directories.extend(dialog.selectedFiles())
            self.file_browser_text.setText(dialog.selectedFiles()[0])
        dialog.deleteLater()

    def browsefilelists(self):
        self.selectfiles("txt")


app = QApplication(sys.argv)
mainwindow = MainWindow()
widget = QtWidgets.QStackedWidget()
widget.addWidget(mainwindow)
widget.setFixedWidth(800)
widget.setFixedHeight(600)
widget.show()
sys.exit(app.exec_())
