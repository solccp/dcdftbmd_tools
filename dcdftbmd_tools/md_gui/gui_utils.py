#!/usr/bin/env python3


import PyQt5.QtWidgets as qtw
from PyQt5.QtCore import QThread, pyqtSignal

class External(QThread):
    finished = pyqtSignal()

    def __init__(self, callback, parent=None):
        super().__init__(parent=parent)
        self.callback = callback
 
    def run(self):
        self.callback()  
        self.finished.emit()

class Actions(qtw.QDialog):
    def __init__(self, parent, label, callback):
        super().__init__(parent)
        self.label = label
        self.initUI()
        self.calc = External(callback)
        self.calc.finished.connect(self.onFinished)
        self.calc.start()
        
    def initUI(self):
        self.setWindowTitle('Progress Bar')
        self.layout = qtw.QVBoxLayout(self)
        self.progress = qtw.QProgressBar()
        
        self.progress.setGeometry(0, 0, 300, 25)
        self.progress.setMinimum(0)
        self.progress.setMaximum(0)
        self.progress.setValue(0)

        self.layout.addWidget(self.progress)
        self.layout.addWidget(qtw.QLabel(self.label))
        self.show()

    def onFinished(self):
        self.done(1)