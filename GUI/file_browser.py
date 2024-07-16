import sys
from PyQt6 import QtWidgets, QtGui
from PyQt6.QtWidgets import QDialog, QApplication, QFileDialog
from PyQt6.uic import loadUi


class FileBrowser(QDialog):
    def __init__(self, parent):
        """
        Initialize the UI from the .ui file, and assign button click values etc

        Args:
            parent (object): Main window that the browser interacts with (to and from)
        """
        super(FileBrowser, self).__init__()
        loadUi("GUI/ui_files/file_dialog.ui", self)
        # locks any other window by making this window modal
        self.setModal(True)
        # disable/enable control buttons
        self.close_filebrowser_window_button.setEnabled(True)
        self.submit_files_button.setEnabled(False)

        # connect buttons to browing results
        self.browse_files_button.clicked.connect(self.browsefiles)
        self.browse_directories_button.clicked.connect(self.browsedirectories)
        self.browse_filelists_button.clicked.connect(self.browsefilelists)

        # connect buttons to deleting selected paths
        self.delete_file_button.clicked.connect(
            lambda: self.delete_item(self.file_list)
        )
        self.delete_directory_button.clicked.connect(
            lambda: self.delete_item(self.directory_list)
        )
        self.delete_filelist_button.clicked.connect(
            lambda: self.delete_item(self.filelist_list)
        )

        # submit files or close window
        self.submit_files_button.clicked.connect(self.submit_files)
        self.close_filebrowser_window_button.clicked.connect(self.close)

        self.parent = parent

    def selectfiles(self, extension: str = ""):
        """
        Method that will allow finding files only with the given extension, and select one or more

        Args:
            extension (str): File extension that is allowed to search for
        """
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
        # enable submit button only if there are files
        if len(self.filelist_list) > 0 or len(self.file_list) > 0:
            self.submit_files_button.setEnabled(True)

    def browsefiles(self):
        """
        Browse results files with file extension given from main UI

        """
        return self.selectfiles(self.parent.filepattern)

    def browsedirectories(self):
        """
        Browses folders and can pick one or more
        """
        dialog = QtWidgets.QFileDialog(
            self,
            options=QFileDialog.Option.DontUseNativeDialog,
            fileMode=QFileDialog.FileMode.Directory,
        )
        # dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
        dialog.setWindowTitle("Choose Directories")

        # dialog.setFileMode()
        for view in dialog.findChildren((QtWidgets.QListView, QtWidgets.QTreeView)):
            if isinstance(view.model(), QtGui.QFileSystemModel):
                view.setSelectionMode(
                    QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
                )

        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.directory_list.insertItems(0, dialog.selectedFiles())
        dialog.deleteLater()

        if len(self.directory_list) > 0:
            self.submit_files_button.setEnabled(True)

    def browsefilelists(self):
        """
        Browses for file lists, which is any file with extension "*.txt"
        """
        return self.selectfiles("*.txt")

    def delete_item(self, list_widget):
        """
        Method to remove items from a QListWidget that are marked in the list

        Args:
            list_widget (QListWidget): List widget from which to take items off
        """
        # get index of the items, and remove from bottom up so index doesn't change
        indices = [x.row() for x in list_widget.selectedIndexes()]
        # then remove each item from bottom up (so indices don't change)
        for index in reversed(indices):
            list_widget.takeItem(index)

    def submit_files(self):
        """
        Submit files to the main window if files are prsent, and updates the text boxes showing how many items were collected in each category
        """
        files = False
        # TODO I can def remove some redundancies here, maybe use lambda expressions or other
        num_dir = str(len(self.directory_list))
        num_fl = str(len(self.file_list))
        num_fll = str(len(self.filelist_list))
        self.parent.numOfDirs.setText(num_dir)
        self.parent.numOfFiles.setText(num_fl)
        self.parent.numOfFilelists.setText(num_fll)

        if len(self.file_list) > 0:
            self.parent.files = self.file_list
            files = True

        if len(self.filelist_list) > 0:
            self.parent.filelists = self.filelist_list
            files = True

        if len(self.directory_list) > 0:
            self.parent.directories = self.directory_list
            files = True

        if files:
            self.parent.submitResultsButton.setEnabled(True)
        else:
            self.parent.submitResultsButton.setEnabled(False)

        self.close_window()

    def close_window(self):
        """
        Close the file browser
        """
        self.close()


if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    filebrowser = FileBrowser()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(filebrowser)
    widget.setFixedWidth(800)
    widget.setFixedHeight(600)
    widget.show()
    sys.exit(app.exec())
