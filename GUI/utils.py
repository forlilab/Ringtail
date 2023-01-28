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

class Interaction:
    def __init__(self):
        self.wanted = False
        self.interaction_type = None
        self.chain = None
        self.res_type = None
        self.res_number = None
        self.atom_name = None
                    
    def set_interaction_type(self, interaction_type):
        if interaction_type != '' and interaction_type != 'Any' and interaction_type != 'None': self.interaction_type = interaction_type
        else: self.interaction_type = None
            
    def set_chain(self, chain):
        if chain != '' and chain != 'Any' and chain != 'None': self.chain = chain
        else: self.chain = None
    
    def set_res_type(self, res_type):
        if res_type != '' and res_type != 'Any' and res_type != 'None': self.res_type = res_type
        else: self.res_type = None
        
    def set_res_number(self, res_number):
        if res_number != '' and res_number != 'Any' and res_number != 'None': self.res_number = res_number
        else: self.res_number = None
        
    def set_atom_name(self, atom_name):
        if atom_name != '' and atom_name != 'Any' and atom_name != 'None': self.atom_name = atom_name
        else: self.atom_name = None
        
    def set_wanted(self, wanted):
        if isinstance(wanted, bool):
            self.wanted = wanted
    
    def __str__(self):
        rep = f"Interaction:{self.interaction_type},Chain:{self.chain},Res_Name:{self.res_type},Res_ID:{self.res_number},Atom_Name:{self.atom_name},Enabled:{str(self.wanted)}"
        return rep


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