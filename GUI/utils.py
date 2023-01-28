from PyQt5 import QtWidgets, QtCore, QtGui

from interaction_widget import Interaction

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

def get_interaction_obj_from_str(s_filter):
    retvalue = Interaction()
    splitted_filter = s_filter.text().split(",")
    retvalue.set_interaction_type(splitted_filter[0].split(':')[-1])
    retvalue.set_chain(splitted_filter[1].split(':')[-1])
    retvalue.set_res_type(splitted_filter[2].split(':')[-1])
    retvalue.set_res_number(splitted_filter[3].split(':')[-1])
    retvalue.set_atom_name(splitted_filter[4].split(':')[-1])
    retvalue.set_wanted(bool(splitted_filter[-1].split(':')[-1]))
    return retvalue

def parse_list_of_items(items):
    retvalue = list()
    for element in items:
        interaction = Interaction()
        interaction.set_interaction_type(element[0].strip())
        interaction.set_chain(element[1].strip())
        interaction.set_res_type(element[2].strip())
        interaction.set_res_number(element[3].strip())
        interaction.set_atom_name(element[4].strip())
        retvalue.append(interaction)
    return retvalue