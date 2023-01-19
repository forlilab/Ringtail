from PyQt5 import QtWidgets, QtCore, QtGui

error_level = [0, 1, 2]

def show_message(message, error_level):
    msg_box = QtWidgets.QMessageBox()
    msg_box.setAutoFillBackground(True)
    
    if error_level == 0:
        msg_box.setIcon(QtWidgets.QMessageBox.Information)
        msg_box.setWindowTitle("Information")
    elif error_level == 1:
        msg_box.setIcon(QtWidgets.QMessageBox.Warning)
        msg_box.setWindowTitle("Warning!")
    elif error_level == 2:
        msg_box.setIcon(QtWidgets.QMessageBox.Critical)
        msg_box.setWindowTitle("Error!")
    
    msg_box.setText(message)
    msg_box.setStandardButtons(QtWidgets.QMessageBox.Ok)
    msg_box.exec_()
    
def browse_file(file_extension=None):
    file = QtWidgets.QFileDialog.getOpenFileName(None, "Find file",
                                                 filter=file_extension)
    return file

def browse_directory():
    dialog = QtWidgets.QDialog()
    directory = QtWidgets.QFileDialog.getExistingDirectory(dialog, "Find Directory",
            QtCore.QDir.currentPath(), )
    return directory

def save_file():
    dialog = QtWidgets.QDialog()
    file = QtWidgets.QFileDialog.getSaveFileName(dialog, "Save file")
    with open(file[0], 'w') as f:
        f.write("File saving try.")
    return file

def get_energy_max_min():
    return (370, -12)

def get_ligands_efficiency_max_min():
    return (3000, -3000)