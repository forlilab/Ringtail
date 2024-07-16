import os
import sys
import util as u
from ringtail_prototype_ui import Ringtail_Prototype_UI

# from PyQt5 import QtWidgets, QtGui
# from PyQt5.QtWidgets import QApplication, QMainWindow
# from PyQt5.uic import loadUi
from PyQt6 import QtWidgets, QtGui
from PyQt6.QtWidgets import QApplication, QMainWindow


class UI_MainWindow(Ringtail_Prototype_UI):
    def __init__(self):
        super(UI_MainWindow, self).__init__()

    def connectUI(self, MainWindow):
        # set up the UI from the main window class
        self.setupUi(MainWindow)


if __name__ == "__main__":
    import sys

    # Create a QT app widget
    app = QtWidgets.QApplication(sys.argv)
    # set the icon that will appear on the dock
    app.setWindowIcon(QtGui.QIcon(":ringtail_head"))
    # create the QT widget main window that will hold and display the app
    MainWindow = QtWidgets.QMainWindow()
    # create a user interface class that will get all the settings
    ui = UI_MainWindow()
    # connect the main window that can be shown, with the ui class that holds all the settings
    ui.connectUI(MainWindow)
    # show the main window
    MainWindow.show()
    # start/execute the application
    sys.exit(app.exec())
