from PyQt6.QtCore import QRunnable, pyqtSlot


# visit www.pyshine.com for more details
class Worker(QRunnable):

    def __init__(self, fnc, *args, **kwargs):
        super(Worker, self).__init__()
        self.fnc = fnc
        self.args = args
        self.kwargs = kwargs

    @pyqtSlot()
    def run(self):

        self.fnc(*self.args, **self.kwargs)
