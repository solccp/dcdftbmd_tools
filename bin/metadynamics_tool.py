#!/usr/bin/env python3

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dcdftbmd_tools.md_gui.models import *
from dcdftbmd_tools.md_gui.viewers import *

import pathlib
import PyQt5.QtWidgets as qtw
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

class MTDTool(qtw.QMainWindow):
    
    def __init__(self):
        super().__init__()
        self.model = MetaDynamicsResultModel()
        self.initUI()
        
        
    def initUI(self):               


        exitAct = qtw.QAction(QIcon('exit.png'), '&Exit', self)        
        exitAct.setShortcut('Ctrl+Q')
        # exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(qtw.qApp.quit)

        


        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        

        #tab
        self.layout = qtw.QVBoxLayout()

        sub_h_layout = qtw.QHBoxLayout()

        sub_h_layout.addWidget(qtw.QLabel('Time:'))
        
        self.cur_time_label = qtw.QLabel('{:8.2f} ps'.format(0.0))
        sub_h_layout.addWidget(self.cur_time_label)
        sub_h_layout.addSpacing(2)
        sub_h_layout.addWidget(qtw.QLabel('Step:'))
        
        self.cur_step_label = qtw.QLabel('0')
        self.total_step_label = qtw.QLabel('0')

        sub_h_layout.addWidget(self.cur_step_label)
        sub_h_layout.addWidget(qtw.QLabel('/'))
        sub_h_layout.addWidget(self.total_step_label)


        self.time_slider = qtw.QSlider(Qt.Horizontal)
        self.time_slider.setEnabled(False)

        sub_h_layout.addWidget(self.time_slider)
        

        main_widget = qtw.QWidget(self)
        main_widget.setLayout(self.layout)

        self.layout.setContentsMargins(1,1,1,1)
        self.layout.addLayout(sub_h_layout)

        self.tab = qtw.QTabWidget()
        
        self.fes_tab = FESWindow()
        self.tab.addTab(self.fes_tab, 'FES')
        
        self.fes_tab_2d = TwoDFESWindow()
        self.tab.addTab(self.fes_tab_2d, '2D FES')

        self.cv_coord_tab = OneDTimeWindow()
        self.tab.addTab(self.cv_coord_tab, 'CV')

        self.cv_height_tab = OneDTimeWindow()
        self.tab.addTab(self.cv_height_tab, 'CV Height')

        self.time_property = TimeCoursePropertyWindow(self.model)
        self.tab.addTab(self.time_property, 'Time vs Property')

        self.layout.addWidget(self.tab)

        self.setCentralWidget(main_widget)


        loadFesAct = qtw.QAction(QIcon(), '&Load Data', self)
        loadFesAct.triggered.connect(self.load_data_dialog)
        fileMenu.addAction(loadFesAct)
        fileMenu.addAction(exitAct)

        # finalizing
        self.setMinimumSize(800, 600)
        self.setWindowTitle('DCDFTBK MetaDynamics Tool')    


        saveAct = qtw.QAction(QIcon.fromTheme('document-save'), 'Save', self)
        saveAct.setShortcut('Ctrl+S')
        saveAct.triggered.connect(self.save_doc)
        
        self.toolbar = self.addToolBar('Save')
        self.toolbar.addAction(saveAct)

        self.show()

    def save_doc(self):      
        fileDialog = qtw.QFileDialog()
        fileDialog.setDirectory(os.path.curdir)
        fileDialog.setFileMode(qtw.QFileDialog.AnyFile)
        fileDialog.setAcceptMode(qtw.QFileDialog.AcceptSave) 
        filename, _ = fileDialog.getSaveFileName()
        if filename != "":
            with open(filename, 'w') as fout:
                if self.tab.currentIndex() == 0:
                    xs, ys = self.fes_tab.curve.getData()
                    print('### FREE ENERGY SURFACE CONSISTING OF  {} GAUSSIANS'.format(int(self.cur_step_label.text())+1), file=fout)
                    for x, y in zip(xs, ys):
                        print('{:20.12f} {:20.12f} {:20.12f}'.format(x, y/627.5095, y),file=fout)


    def load_data(self, folder):
        self.folder = pathlib.Path(folder)
        self.model.set_base_folder(self.folder)

        if (len(self.model.get_gau_pots()) > 0):
            self.time_slider.setEnabled(True)
            self.time_slider.valueChanged.connect(self.update_plots)

            last_step = (len(self.model.gau_pots)-1)

            self.cur_step_label.setText(str(last_step))
            self.total_step_label.setText(str(last_step))

            self.time_slider.setMaximum(last_step)
            self.time_slider.setValue(last_step)

            cv_coord_model = CVCoordAdpater(self.model)
            self.cv_coord_tab.set_model(cv_coord_model)

            cv_height_model = CVHeightAdpater(self.model)
            self.cv_height_tab.set_model(cv_height_model)

            self.init_plots()
            self.update_plots(last_step)

            if self.model.get_fes_dimension() == 1:
                self.tab.setCurrentIndex(0)
                self.tab.setTabEnabled(0, True)
                self.tab.setTabEnabled(1, False)
            elif self.model.get_fes_dimension() == 2:
                self.tab.setCurrentIndex(1)
                self.tab.setTabEnabled(0, False)
                self.tab.setTabEnabled(1, True)

            
    def clean_data(self):
        self.model = MetaDynamicsResultModel()
        self.time_slider.setEnabled(False)

    def load_data_dialog(self):
        dialog = qtw.QFileDialog(self, 'Load Data Folder', pathlib.os.curdir)
        folder = dialog.getExistingDirectory()
        
        if folder != "":
            self.load_data(folder)
            
        
    def init_plots(self):
        self.time_property.update_plot()

    def update_plots(self, step):
        self.cur_step_label.setText(str(step))
        self.cur_time_label.setText('{:8.2f} ps'.format((step+1)*self.model.gaussian_interval_time/1000.0))
        
        if self.model.get_fes_dimension() == 1:
            self.fes_tab.plot(step, self.model)
            
        if self.model.get_fes_dimension() == 2:
            self.fes_tab_2d.plot(step, self.model)
        
        self.cv_coord_tab.update_time(step)
        self.cv_height_tab.update_time(step)
        self.time_property.update_time(step)
        
        

if __name__ == '__main__':
    
    app = qtw.QApplication(sys.argv)
    tool = MTDTool()
    tool.resize(1024, 768)
    # try:
    tool.load_data(os.curdir)
    # except Exception as e:
        # print(e)
        # tool.clean_data()
    sys.exit(app.exec_())

