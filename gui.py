
import sys
from os import system
from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit,
    QTextEdit, QGridLayout, QApplication, QPushButton, )


class Example(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()


    def initUI(self):

        cnv_cases = QLabel('Enter the number of CNVs observed in cases:')
        tot_cases = QLabel('Enter the total number of patients in cases:')
        cnv_ctrl = QLabel('Enter the number of CNVs observed in controls:')
        tot_ctrl = QLabel('Enter the total number of patients in controls:')

        cnv_casesEdit = QLineEdit()
        tot_casesEdit = QLineEdit()
        cnv_ctrlEdit = QLineEdit()
        tot_ctrlEdit = QLineEdit()
        btnEdit = QPushButton('Calculate Penetrance')

        grid = QGridLayout()
        grid.setSpacing(10)

        grid.addWidget(cnv_cases, 1, 0)
        grid.addWidget(cnv_casesEdit, 1, 1)

        grid.addWidget(tot_cases, 2, 0)
        grid.addWidget(tot_casesEdit, 2, 1)

        grid.addWidget(cnv_ctrl, 3, 0)
        grid.addWidget(cnv_ctrlEdit, 3, 1)

        grid.addWidget(tot_ctrl, 4, 0)
        grid.addWidget(tot_ctrlEdit, 4, 1)

        grid.addWidget(btnEdit,5,0)

        btnEdit.clicked.connect(self.doAction)

        self.setLayout(grid)

        self.setGeometry(500, 500, 500, 200)
        self.setWindowTitle('Penetrance')
        self.show()

    def doAction(self):
        system("python v1.py");

if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
