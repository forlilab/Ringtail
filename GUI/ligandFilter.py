from PyQt5 import QtCore, QtGui, QtWidgets
from utils import LigandFilter, show_message


class Ui_Dialog(QtWidgets.QDialog):
    def __init__(self, wanted=None):
        super(Ui_Dialog, self).__init__()
        self.substructure = False
        self.ligand_name = False
        self.ligand_filter = LigandFilter()
        self.ligand_filter.enabled = wanted
        
        self.setupUi(self)
    
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(450, 287)
        self.gridLayout = QtWidgets.QGridLayout(Dialog)
        self.gridLayout.setObjectName("gridLayout")
        self.ligandNameRadioButton = QtWidgets.QRadioButton(Dialog)
        self.ligandNameRadioButton.setObjectName("ligandNameRadioButton")
        self.gridLayout.addWidget(self.ligandNameRadioButton, 0, 0, 1, 3)
        self.ligandSubstructureRadioButton = QtWidgets.QRadioButton(Dialog)
        self.ligandSubstructureRadioButton.setObjectName("ligandSubstructureRadioButton")
        self.gridLayout.addWidget(self.ligandSubstructureRadioButton, 0, 5, 1, 2)
        self.ligandNameLabel = QtWidgets.QLabel(Dialog)
        self.ligandNameLabel.setObjectName("ligandNameLabel")
        self.gridLayout.addWidget(self.ligandNameLabel, 1, 0, 1, 1)
        self.ligandNameLineEdit = QtWidgets.QLineEdit(Dialog)
        self.ligandNameLineEdit.setObjectName("ligandNameLineEdit")
        self.gridLayout.addWidget(self.ligandNameLineEdit, 1, 1, 1, 5)
        self.substructureLabel = QtWidgets.QLabel(Dialog)
        self.substructureLabel.setObjectName("substructureLabel")
        self.gridLayout.addWidget(self.substructureLabel, 2, 0, 1, 2)
        self.substructureLineEdit = QtWidgets.QLineEdit(Dialog)
        self.substructureLineEdit.setObjectName("substructureLineEdit")
        self.gridLayout.addWidget(self.substructureLineEdit, 2, 3, 1, 3)
        self.coordinatesCheckBox = QtWidgets.QCheckBox(Dialog)
        self.coordinatesCheckBox.setObjectName("coordinatesCheckBox")
        self.gridLayout.addWidget(self.coordinatesCheckBox, 3, 0, 1, 4)
        self.coordinatesLabel = QtWidgets.QLabel(Dialog)
        self.coordinatesLabel.setObjectName("coordinatesLabel")
        self.gridLayout.addWidget(self.coordinatesLabel, 4, 0, 1, 3)
        self.cutoffLabel = QtWidgets.QLabel(Dialog)
        self.cutoffLabel.setObjectName("cutoffLabel")
        self.gridLayout.addWidget(self.cutoffLabel, 4, 5, 1, 1)
        self.indexLabel = QtWidgets.QLabel(Dialog)
        self.indexLabel.setObjectName("indexLabel")
        self.gridLayout.addWidget(self.indexLabel, 4, 6, 1, 1)
        self.xDoubleSpinBox = QtWidgets.QDoubleSpinBox(Dialog)
        self.xDoubleSpinBox.setObjectName("xDoubleSpinBox")
        self.gridLayout.addWidget(self.xDoubleSpinBox, 5, 0, 1, 2)
        self.yDoubleSpinBox = QtWidgets.QDoubleSpinBox(Dialog)
        self.yDoubleSpinBox.setObjectName("yDoubleSpinBox")
        self.gridLayout.addWidget(self.yDoubleSpinBox, 5, 2, 1, 2)
        self.zDoubleSpinBox = QtWidgets.QDoubleSpinBox(Dialog)
        self.zDoubleSpinBox.setObjectName("zDoubleSpinBox")
        self.gridLayout.addWidget(self.zDoubleSpinBox, 5, 4, 1, 1)
        self.cutoffDoubleSpinBox = QtWidgets.QDoubleSpinBox(Dialog)
        self.cutoffDoubleSpinBox.setObjectName("cutoffDoubleSpinBox")
        self.gridLayout.addWidget(self.cutoffDoubleSpinBox, 5, 5, 1, 1)
        self.indexSpinBox = QtWidgets.QSpinBox(Dialog)
        self.indexSpinBox.setObjectName("indexSpinBox")
        self.gridLayout.addWidget(self.indexSpinBox, 5, 6, 1, 1)
        self.okCancelButtonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.okCancelButtonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.okCancelButtonBox.setObjectName("okCancelButtonBox")
        self.gridLayout.addWidget(self.okCancelButtonBox, 6, 5, 1, 2)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
        
        self.name_elements = [self.ligandNameLabel,
                              self.ligandNameLineEdit]
        
        self.sub_elements = [self.substructureLabel,
                             self.substructureLineEdit,
                             self.coordinatesCheckBox]
        
        self.coordinates_elements = [self.coordinatesLabel,
                                     self.cutoffLabel,
                                     self.indexLabel,
                                     self.xDoubleSpinBox,
                                     self.yDoubleSpinBox,
                                     self.zDoubleSpinBox,
                                     self.cutoffDoubleSpinBox,
                                     self.indexSpinBox]
        
        # EVENTS
        self.ligandNameRadioButton.clicked.connect(self.set_mode)
        self.ligandSubstructureRadioButton.clicked.connect(self.set_mode)
        
        # Name
        self.ligandNameLineEdit.textChanged.connect(self.set_ligand_name)
        
        # Substructure
        self.substructureLineEdit.textChanged.connect(self.set_substructure_match)
        self.coordinatesCheckBox.stateChanged.connect(self.enable_coordinates)
        self.xDoubleSpinBox.valueChanged.connect(self.set_x_coord)
        self.yDoubleSpinBox.valueChanged.connect(self.set_y_coord)
        self.zDoubleSpinBox.valueChanged.connect(self.set_z_coord)
        self.cutoffDoubleSpinBox.valueChanged.connect(self.set_cutoff)
        self.indexSpinBox.valueChanged.connect(self.set_index)
        
        self.okCancelButtonBox.accepted.connect(self.get_ligand_filter_obj)
        self.okCancelButtonBox.rejected.connect(self.reject)
        
        self.ligandNameRadioButton.click()
        self.xDoubleSpinBox.clear()
        self.yDoubleSpinBox.clear()
        self.zDoubleSpinBox.clear()
        self.cutoffDoubleSpinBox.clear()
        self.indexSpinBox.clear()
    
    def set_mode(self):
        self.okCancelButtonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(False)
        if self.ligandNameRadioButton.isChecked():
            self.substructure = False
            self.ligand_name = True
            self.enable_substructure_elements(False)
            self.ligandNameLineEdit.clear()
            self.substructureLineEdit.clear()
            self.xDoubleSpinBox.clear()
            self.yDoubleSpinBox.clear()
            self.zDoubleSpinBox.clear()
            self.cutoffDoubleSpinBox.clear()
            self.indexSpinBox.clear()
        else:
            self.substructure = True
            self.ligand_name = False
            self.enable_substructure_elements(True)
            self.ligandNameLineEdit.clear()
            self.substructureLineEdit.clear()
            self.xDoubleSpinBox.clear()
            self.yDoubleSpinBox.clear()
            self.zDoubleSpinBox.clear()
            self.cutoffDoubleSpinBox.clear()
            self.indexSpinBox.clear()
            self.enable_coordinates()
            self.coordinatesCheckBox.setEnabled(False)
            
    def enable_coordinates(self):
        if self.coordinatesCheckBox.isChecked():
            self.ligand_filter.include_coordinates = True
            for element in self.coordinates_elements:
                element.setEnabled(True)
        else:
            self.ligand_filter.include_coordinates = False
            for element in self.coordinates_elements:
                element.setEnabled(False)
    
    def enable_substructure_elements(self, enabled):
        for element in self.name_elements:
            element.setVisible(not enabled)
            element.setEnabled(not enabled)
        for element in self.sub_elements:
            element.setVisible(enabled)
            element.setEnabled(enabled)
        for element in self.coordinates_elements:
            element.setVisible(enabled)
                
    def set_ligand_name(self):
        if self.ligand_name:
            if self.ligandNameLineEdit.text() != '':
                self.ligand_filter.set_ligand_name(self.ligandNameLineEdit.text().strip())
                self.okCancelButtonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(True)
            else:
                self.okCancelButtonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(False)
                
    def set_substructure_match(self):
        if self.substructure:
            if self.substructureLineEdit.text() != '':
                self.okCancelButtonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(True)
                self.ligand_filter.set_substructure_match(self.substructureLineEdit.text())
                self.coordinatesCheckBox.setEnabled(True)
            else:
                self.okCancelButtonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(False)
                self.coordinatesCheckBox.setChecked(False)
                self.coordinatesCheckBox.setEnabled(False)
    
    def set_x_coord(self):
        if self.substructure:
            self.ligand_filter.set_x(self.xDoubleSpinBox.value())
    
    def set_y_coord(self):
        if self.substructure:
            self.ligand_filter.set_y(self.yDoubleSpinBox.value())
    
    def set_z_coord(self):
        if self.substructure:
            self.ligand_filter.set_z(self.zDoubleSpinBox.value())
    
    def set_cutoff(self):
        if self.substructure:
            self.ligand_filter.set_cutoff(self.cutoffDoubleSpinBox.value())
    
    def set_index(self):
        if self.substructure:
            self.ligand_filter.set_index(self.indexSpinBox.value())
            
    def get_ligand_filter_obj(self):
        if self.ligand_name:
            self.ligand_filter.set_ligand_name(self.ligandNameLineEdit.text())
            self.ligand_filter.set_substructure_match(None)
            self.ligand_filter.set_coordinates(None, None, None)
            self.ligand_filter.set_cutoff(None)
            self.ligand_filter.set_index(None)
            print(self.ligand_filter)
            self.accept()
        else:
            self.ligand_filter.set_ligand_name(None)
            self.ligand_filter.set_substructure_match(self.substructureLineEdit.text())
            if self.coordinatesCheckBox.isChecked():
                if self.coordinates_sanity_check():
                    print(self.ligand_filter)
                    self.accept()
                else: show_message("Error if coordinates checkbox is active all the fields must be filled in.",
                                   2)
            else:
                self.ligand_filter.set_coordinates(None, None, None)
                self.ligand_filter.set_cutoff(None)
                self.ligand_filter.set_index(None)
                print(self.ligand_filter)
                self.accept()

    def coordinates_sanity_check(self):
        return self.xDoubleSpinBox.value() != 0.0 and self.xDoubleSpinBox.text() != '' \
            and self.yDoubleSpinBox.value() != 0.0 and self.yDoubleSpinBox.text() != '' \
            and self.zDoubleSpinBox.value() != 0.0  and self.zDoubleSpinBox.text() != '' \
            and self.cutoffDoubleSpinBox.value() > 0.0 and self.indexSpinBox.text() != ''
    
    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.ligandNameRadioButton.setText(_translate("Dialog", "Ligand Name"))
        self.ligandSubstructureRadioButton.setText(_translate("Dialog", "Substructure Match"))
        self.ligandNameLabel.setText(_translate("Dialog", "Name:"))
        self.substructureLabel.setText(_translate("Dialog", "Substructure:"))
        self.coordinatesCheckBox.setText(_translate("Dialog", "Include coordinates"))
        self.coordinatesLabel.setText(_translate("Dialog", "Coordinates (x, y, z):"))
        self.cutoffLabel.setText(_translate("Dialog", "Cutoff:"))
        self.indexLabel.setText(_translate("Dialog", "Index:"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ui = Ui_Dialog()
    ui.show()
    sys.exit(app.exec_())
