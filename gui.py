import sys
from pyQt4 import QtGui

app = QtGui.QApplcation(sys.argv)

class Window(QtGui.QMainWindow):

    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("Penetrance")
        # self.setWindowIcon(QtGui.QIcon('logo.png'))
        self.show()

GUI = Window()
sys.exit(app.exec_())
