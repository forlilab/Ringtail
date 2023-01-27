# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'UI/interaction.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
from list_example import list_example

class Interaction:
    def __init__(self):
        self.wanted = False
        self.interaction_type = None
        self.chain = None
        self.res_type = None
        self.res_number = None
        self.atom_name = None
        
    def set_interaction_type(self, interaction_type):
        self.interaction_type = interaction_type
    
    def set_chain(self, chain):
        self.chain = chain
    
    def set_res_type(self, res_type):
        self.res_type = res_type
        
    def set_res_number(self, res_number):
        self.res_number = res_number
        
    def set_atom_name(self, atom_name):
        self.atom_name = atom_name
    
    def set_wanted(self, wanted):
        if isinstance(wanted, bool):
            self.wanted = wanted
    
    def __str__(self):
        rep = f"Interaction: {self.interaction_type}\nChain: {self.chain}\nRes Name: {self.res_type}\nRes ID: {self.res_number}\nAtom Name: {self.atom_name}"
        return rep
    
class Ui_Form(QtWidgets.QWidget):
    def __init__(self, list_of_available_items, parent : QtWidgets.QWidget):
        super().__init__(parent)
        self.interactions_type = ['R', 'H', 'vdW']
        self.items = parse_list_of_items(list_of_available_items)
        self.filtered_items = self.items
        self.filter = Interaction()
        self.setupUi(parent)
        
        self.interaction_filter = False
        self.chain_filter = False
        self.res_name_filter = False
        self.res_id_filter = False
        self.atom_name_filter = False
    
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(542, 128)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.gridLayout = QtWidgets.QGridLayout(Form)
        self.gridLayout.setObjectName("gridLayout")
        self.interactionTypeLabel = QtWidgets.QLabel(Form)
        self.interactionTypeLabel.setFocusPolicy(QtCore.Qt.NoFocus)
        self.interactionTypeLabel.setObjectName("interactionTypeLabel")
        self.gridLayout.addWidget(self.interactionTypeLabel, 0, 0, 1, 2)
        self.interactionTypeComboBox = QtWidgets.QComboBox(Form)
        self.interactionTypeComboBox.setObjectName("interactionTypeComboBox")
        self.gridLayout.addWidget(self.interactionTypeComboBox, 0, 2, 1, 2)
        self.wantedRadioButton = QtWidgets.QRadioButton(Form)
        self.wantedRadioButton.setObjectName("wantedRadioButton")
        self.gridLayout.addWidget(self.wantedRadioButton, 0, 4, 1, 1)
        self.unwantedRadioButton = QtWidgets.QRadioButton(Form)
        self.unwantedRadioButton.setObjectName("unwantedRadioButton")
        self.gridLayout.addWidget(self.unwantedRadioButton, 0, 5, 1, 1)
        self.interactionChainLabel = QtWidgets.QLabel(Form)
        self.interactionChainLabel.setObjectName("interactionChainLabel")
        self.gridLayout.addWidget(self.interactionChainLabel, 1, 0, 1, 1)
        self.interactionResTypeLabel = QtWidgets.QLabel(Form)
        self.interactionResTypeLabel.setObjectName("interactionResTypeLabel")
        self.gridLayout.addWidget(self.interactionResTypeLabel, 1, 1, 1, 2)
        self.interactionResNLabel = QtWidgets.QLabel(Form)
        self.interactionResNLabel.setObjectName("interactionResNLabel")
        self.gridLayout.addWidget(self.interactionResNLabel, 1, 3, 1, 1)
        self.interactionAtomNameLabel = QtWidgets.QLabel(Form)
        self.interactionAtomNameLabel.setObjectName("interactionAtomNameLabel")
        self.gridLayout.addWidget(self.interactionAtomNameLabel, 1, 4, 1, 2)
        self.interactionChainComboBox = QtWidgets.QComboBox(Form)
        self.interactionChainComboBox.setObjectName("interactionChainComboBox")
        self.gridLayout.addWidget(self.interactionChainComboBox, 2, 0, 1, 1)
        self.interactionResTypeComboBox = QtWidgets.QComboBox(Form)
        self.interactionResTypeComboBox.setObjectName("interactionResTypeComboBox")
        self.gridLayout.addWidget(self.interactionResTypeComboBox, 2, 1, 1, 2)
        self.interactionResNComboBox = QtWidgets.QComboBox(Form)
        self.interactionResNComboBox.setObjectName("interactionResNComboBox")
        self.gridLayout.addWidget(self.interactionResNComboBox, 2, 3, 1, 1)
        self.interactionAtomNameComboBox = QtWidgets.QComboBox(Form)
        self.interactionAtomNameComboBox.setObjectName("interactionAtomNameComboBox")
        self.gridLayout.addWidget(self.interactionAtomNameComboBox, 2, 4, 1, 2)
        self.interactionDialogButtonBox = QtWidgets.QDialogButtonBox(Form)
        self.interactionDialogButtonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.interactionDialogButtonBox.setObjectName("interactionDialogButtonBox")
        self.gridLayout.addWidget(self.interactionDialogButtonBox, 3, 4, 1, 2)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

        
        self.interactionTypeComboBox.currentIndexChanged.connect(self.set_interaction_type)
        self.wantedRadioButton.clicked.connect(self.set_filter_wanted)
        self.unwantedRadioButton.clicked.connect(self.set_filter_unwanted)
        
        self.interactionTypeComboBox.currentIndexChanged.connect(self.update_from_interaction_type)
        self.interactionChainComboBox.currentIndexChanged.connect(self.update_from_chain)
        # self.interactionResTypeComboBox.currentIndexChanged.connect(None)
        # self.interactionResNComboBox.currentIndexChanged.connect(None)
        # self.interactionAtomNameComboBox.currentIndexChanged.connect(None)
        
        # DEFAULTS
        self.populate_comboboxes(self.items)
        self.wantedRadioButton.click()
        
    # CLASS METHODS
    def populate_comboboxes(self, 
                            items, 
                            interaction=False, 
                            chain=False, 
                            res_type=False, 
                            res_number=False, 
                            atom_name=False):
        
        if not interaction:
            self.interactionTypeComboBox.clear()
            self.interactionTypeComboBox.addItems(str(t) for t in set([x.interaction_type for x in items]))
        if not chain:
            self.interactionChainComboBox.clear()
            self.interactionChainComboBox.addItems(str(t) for t in set([x.chain for x in items]))
        if not res_type:
            self.interactionResTypeComboBox.clear()
            self.interactionResTypeComboBox.addItems(str(t) for t in set([x.res_type for x in items]))
        if not res_number:
            self.interactionResNComboBox.clear()
            self.interactionResNComboBox.addItems(str(t) for t in set([x.res_number for x in items]))
        if not atom_name:
            self.interactionAtomNameComboBox.clear()
            self.interactionAtomNameComboBox.addItems(str(t) for t in set([x.atom_name for x in items]))
    
    def update_from_interaction_type(self):
        # The interaction filter is the first filter applied
        self.interaction_filter = True
        self.chain_filter = False
        self.res_name_filter = False
        self.res_id_filter = False
        self.atom_name_filter = False
        
        interaction_type = self.interactionTypeComboBox.currentText()
        self.filtered_items = list(filter(lambda x: x.interaction_type == interaction_type, self.items))
        self.populate_comboboxes(self.filtered_items, interaction=self.interaction_filter)
    
    def set_interaction_type(self):
        self.filter.set_interaction_type(self.interactionTypeComboBox.currentText())
        
    def set_filter_wanted(self):
        if self.wantedRadioButton.isChecked():
            self.filter.set_wanted(True)
        else:
            self.filter.set_wanted(False)
    
    def set_filter_unwanted(self):
        if self.unwantedRadioButton.isChecked():
            self.filter.set_wanted(False)
        else:
            self.filter.set_wanted(True)
    
    def get_interaction_obj(self):
        return self.filter

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.interactionTypeLabel.setText(_translate("Form", "Select interaction type:"))
        self.wantedRadioButton.setText(_translate("Form", "Enable"))
        self.unwantedRadioButton.setText(_translate("Form", "Disable"))
        self.interactionChainLabel.setText(_translate("Form", "Chain:"))
        self.interactionResTypeLabel.setText(_translate("Form", "Res Type:"))
        self.interactionResNLabel.setText(_translate("Form", "Res #:"))
        self.interactionAtomNameLabel.setText(_translate("Form", "Atom Name:"))

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


if __name__ == "__main__":
    interactions = parse_list_of_items(list_example)
    print(f"Unique Interactions: {set([x.interaction_type for x in interactions])}")
    print(f"Unique Chains: {set([x.chain for x in interactions])}")
    print(f"Unique Res Names: {set([x.res_type for x in interactions])}")
    print(f"Unique Res IDs: {set([x.res_number for x in interactions])}")
    print(f"Unique Atom Names: {set([x.atom_name for x in interactions])}")