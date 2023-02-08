from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
# CLASSES
class MyDelegate(QItemDelegate):
    def __init__(self, parent=None, *args):
        QItemDelegate.__init__(self, parent, *args)

    def paint(self, painter, option, index):
        painter.save()

        # set background color
        painter.setPen(QPen(Qt.NoPen))
        if option.state & QStyle.State_Selected:
            # If the item is selected, always draw background red
            painter.setBrush(QBrush(Qt.transparent))
        else:
            c = index.data(Qt.DisplayRole+1) # Get the color
            painter.setBrush(QBrush(QColor(c)))

        # Draw the background rectangle            
        painter.drawRect(option.rect)

        # Draw the bottom border
        # option.rect is the shape of the item; top left bottom right
        # e.g. 0, 0, 256, 16 in the parent listwidget
        painter.setPen(QPen(Qt.black))        
        painter.drawLine(option.rect.bottomLeft(), option.rect.bottomRight())

        # Draw the text
        painter.setPen(QPen(Qt.black))
        text = index.data(Qt.DisplayRole)
        # Adjust the rect (to pad)
        option.rect.setLeft(5)
        option.rect.setRight(option.rect.right()-5)
        painter.drawText(option.rect, Qt.AlignLeft, text)

        painter.restore()

class Interaction:
    def __init__(self):
        self.wanted = False
        self.interaction_type = None
        self.chain = None
        self.res_type = None
        self.res_number = None
        self.atom_name = None
        self.accepted = False
        self.enabled = True
                    
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
    
    def set_enabled(self, enabled):
        if isinstance(enabled, bool):
            self.enabled = enabled
    
    def __str__(self):
        rep = f"Interaction:{self.interaction_type},Chain:{self.chain},Res_Name:{self.res_type},Res_ID:{self.res_number},Atom_Name:{self.atom_name},Wanted:{str(self.wanted)}"
        return rep

class LigandFilter:
    def __init__(self):
        self.ligand_name = None
        self.substructure_match = None
        self.include_coordinates = False
        self.x = None
        self.y = None
        self.z = None
        self.cutoff = None
        self.index = None
        self.enabled = True
        
    def set_ligand_name(self, name):
        self.ligand_name = name
    
    def set_substructure_match(self, sub):
        self.substructure_match = sub
    
    def set_coordinates(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
            
    def set_x(self, x):
        self.x = x
    
    def set_y(self, y):
        self.y = y
    
    def set_z(self, z):
        self.z = z
            
    def set_cutoff(self, cutoff):
        self.cutoff = cutoff
    
    def set_index(self, idx):
        self.index = idx
            
    def set_wanted(self, enabled):
        if isinstance(enabled, bool):
            self.enabled = enabled
            
    def __str__(self):
        if self.ligand_name is not None:
            rep = f"Name: {self.ligand_name}"
            rep += f"\nWanted: {self.enabled}"
            return rep
        elif self.substructure_match is not None:
            rep = f"Sub: {self.substructure_match}"
            rep += f"\nWanted: {self.enabled}"
            if self.include_coordinates:
                rep += f"\nx:{self.x}, y:{self.y}, z:{self.x}, cutoff:{self.cutoff}, idx:{self.index}"
            return rep

error_level = [0, 1, 2]

def show_message(message, error_level):
    msg_box = QMessageBox()
    msg_box.setAutoFillBackground(True)
    
    if error_level == 0:
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setWindowTitle("Information")
    elif error_level == 1:
        msg_box.setIcon(QMessageBox.Warning)
        msg_box.setWindowTitle("Warning!")
    elif error_level == 2:
        msg_box.setIcon(QMessageBox.Critical)
        msg_box.setWindowTitle("Error!")
    
    msg_box.setText(message)
    msg_box.setStandardButtons(QMessageBox.Ok)
    msg_box.exec_()
    
def browse_file(file_extension=None):
    file = QFileDialog.getOpenFileName(None, "Find file",
                                                 filter=file_extension)
    return file

def browse_directory():
    dialog = QDialog()
    directory = QFileDialog.getExistingDirectory(dialog, "Find Directory",
            QDir.currentPath(), )
    return directory

def save_file():
    dialog = QDialog()
    file = QFileDialog.getSaveFileName(dialog, "Save file")
    with open(file[0], 'w') as f:
        f.write("File saving try.")
    return file

def get_energy_max_min():
    return (370, -12)

def get_ligands_efficiency_max_min():
    return (3000, -3000)

def get_ligand_obj_from_str(s_filter):
    retvalue = LigandFilter()
    set_name = False
    splitted_filter = s_filter.text().split("\n")
    retvalue.set_wanted(bool(splitted_filter[1].split('Wanted:')[-1].strip()))
    for x in splitted_filter:
        if x.find('Name:') >= 0:
            set_name = True
    if set_name is True:
        retvalue.set_ligand_name(splitted_filter[0].split('Name:')[-1].strip())
    else:
        retvalue.set_substructure_match(splitted_filter[0].split('Sub:')[-1].strip())
    if len(splitted_filter) > 2:
        retvalue.include_coordinates = True
        splitted_coordinates = splitted_filter[-1].split(',')
        x = float(splitted_coordinates[0].split('x:')[-1].strip())
        y = float(splitted_coordinates[1].split('y:')[-1].strip())
        z = float(splitted_coordinates[2].split('z:')[-1].strip())
        cutoff = float(splitted_coordinates[3].split('cutoff:')[-1].strip())
        idx = int(splitted_coordinates[-1].split('idx:')[-1].strip())
        retvalue.set_coordinates(x, y, z)
        retvalue.set_cutoff(cutoff)
        retvalue.set_index(idx)
    return retvalue


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
            
    