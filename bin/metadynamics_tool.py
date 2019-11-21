#!/usr/bin/env python3

import sys
import PyQt5.QtWidgets as qtw
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt
import numpy as np
import matplotlib
import matplotlib.pylab as plt
from scipy.signal import find_peaks

import pathlib
matplotlib.use('QT5Agg')
matplotlib.rcParams["figure.dpi"] = 200.0
matplotlib.rcParams["figure.figsize"] = [10, 6]

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


class CVCoordWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.model = None
        self.line = None

    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        self.figure = plt.figure()
        self.figure.subplots_adjust(top=0.950, bottom=0.15, left=0.1, right=0.95)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.axis = self.figure.add_subplot(1,1,1)
        # self.figure.tight_layout()
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)

    def plot(self, step, model):
        if (self.model == model):
            pass
        else:
            self.model = model

            self.axis.cla()
            
            xs = [self.model.gaussian_interval_time*x/1000.0 for x in range(len(self.model.gau_pots))]
            for i in range(len(self.model.gau_pots[0])):
                values = [x[i].cv_coord for x in self.model.gau_pots]
                self.axis.plot(xs, values, label='CV{}'.format(i+1))
            
            self.axis.set_xlabel('Time [ps]')
            self.axis.set_ylabel('CV1 [{}]'.format(self.model.cvs[0]))
            self.ylim = self.axis.get_ylim()
            self.axis.legend()
             
        if (self.line):
            self.axis.lines.pop()
            self.line = None

        x_pos = self.model.gaussian_interval_time*(step+1)/1000.0
        
        
        self.line = self.axis.plot([x_pos, x_pos], [self.ylim[0], self.ylim[1]], '--', color='red')
        self.canvas.draw()


class CVHeightWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.model = None
        self.line = None

    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.figure.subplots_adjust(top=0.950, bottom=0.15, left=0.1, right=0.95)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.axis = self.figure.add_subplot(1,1,1)
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)

    def plot(self, step, model):
        if (self.model == model):
            pass
        else:

            self.model = model

            self.axis.cla()
            
            xs = [self.model.gaussian_interval_time*x/1000.0 for x in range(len(self.model.gau_pots))]
            for i in range(len(self.model.gau_pots[0])):
                values = [x[i].height*627.5095 for x in self.model.gau_pots]
                self.axis.plot(xs, values, label='CV{}'.format(i+1))
            
            self.axis.set_xlabel('Time [ps]')
            self.axis.set_ylabel('Gaussian Height of CV1 [kcal/mol]')
            self.ylim = self.axis.get_ylim()
            self.axis.legend()
            
             
        if (self.line):
            self.axis.lines.pop()
            self.line = None
    
        x_pos = self.model.gaussian_interval_time*(step+1)/1000.0
        
        self.line = self.axis.plot([x_pos, x_pos], [self.ylim[0], self.ylim[1]], '--', color='red')
        self.canvas.draw()

class FESWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.axis = self.figure.add_subplot(1,1,1)
        self.figure.subplots_adjust(top=0.950, bottom=0.15, left=0.1, right=0.95)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)

    def find_minima_barriers(self, coords, values):
        maxima, _ = find_peaks(values)
        minima, _ = find_peaks([-x for x in values])

        return (maxima, minima)

    def plot(self, step, model):
        coords = model.fes[step].coords
        values = [x*627.5095 for x in model.fes[step].values]
      
        maxima, minima = self.find_minima_barriers(coords, values)

        
        # if 1d
        if (len(coords[0]) == 1):
            xs = [x[0] for x in coords]
            
            self.axis.cla()
            self.axis.plot(xs, values)

            self.axis.set_xlabel(model.cvs[0])
            self.axis.set_ylabel('Free energy [kcal/mol]')
  

            index = zip(minima,minima[1:])
            for step, ind in enumerate(index):
                ind_min_1, ind_min_2 = ind
                ind_max = maxima[step]
             
                x_min = xs[ind_min_1]
                x_max = xs[ind_min_2]
                x_width = x_max-x_min
                y_max = values[ind_max]
                self.axis.plot([x_min, x_max], [y_max, y_max], '--' )

                min_1_x_pos = xs[ind_min_1]
                min_1_y_pos = values[ind_min_1]
                min_2_x_pos = xs[ind_min_2]
                min_2_y_pos = values[ind_min_2]

                energy_1 = y_max - min_1_y_pos
                energy_2 = y_max - min_2_y_pos
    
                self.axis.annotate('', xy=(min_1_x_pos, min_1_y_pos), xytext=(min_1_x_pos, y_max), arrowprops=dict(arrowstyle='<|-|>'))
                
                self.axis.text(min_1_x_pos+x_width*0.01, (min_1_y_pos+y_max)/2.0, '{:.2f}'.format(energy_1), fontsize=9, 
                horizontalalignment='left', verticalalignment='center', bbox=dict(boxstyle='square,pad=0.1', fc='w', ec='none'))

                self.axis.annotate('', xy=(min_2_x_pos, min_2_y_pos), xytext=(min_2_x_pos, y_max), arrowprops=dict(arrowstyle='<|-|>'))
                self.axis.text(min_2_x_pos-x_width*0.01, (min_2_y_pos+y_max)/2.0, '{:.2f}'.format(energy_2), fontsize=9, 
                horizontalalignment='right', verticalalignment='center', bbox=dict(boxstyle='square,pad=0.1', fc='w', ec='none'))           
            
            
            self.canvas.draw()
        else:
            raise RuntimeError('More than 1D FES plotting is not implemented')


class TimeCoursePropertyWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()

    

class GaussianPotentialModel:
    def __init__(self, height, width, cv_num, cv_coord):
        super().__init__()            
        self.height = height
        self.width = width
        self.cv_num = cv_num
        self.cv_coord = cv_coord

class FESPotential:
    def __init__(self):
        super().__init__()
        self.coords = []
        self.values = []
        


class MetaDynamicsResultModel:
    
    def __init__(self):
        super().__init__()
        self.gau_pots = []
        self.fes = []
        self.cvs = []
    
    def parseFolder(self, folderName):
        self.loadData_metacvs(folderName.joinpath('metacv.dat'))
        self.loadData_fes(folderName.joinpath('fes.dat'))
        self.loadData_biaspot(folderName.joinpath('biaspot'))

    def loadData_metacvs(self, fileName):
        with open(fileName, 'r') as f:
            line = next(f)
            if ('RATIONALCOORDINATIONNUMBER' in line):
                self.cvs.append('Coordination Numbers')



    """ GAUSSIAN BIAS POTENTIAL:        1
    *** AT T=         80.00 FSEC, THIS RUN'S STEP NO.=      80
     Gaussian height     =        0.0010000000 a.u.
     Collective variable =                   1
       Coordinate        =        0.0208250915
       Gaussian width    =        0.1000000000
    """
    def loadData_biaspot(self, fileName):
        with open(fileName, 'r') as f:
            line = next(f)
            
            while('GAUSSIAN BIAS POTENTIAL:' in line):
                local_gau_pots = []
                arr = line.split()
                pot_step = int(arr[3])
                title_line = next(f)
                arr = title_line.split()
                interval_step = int(arr[9])
                if (pot_step == 1):
                    self.gaussian_interval_step = interval_step
                    self.gaussian_interval_time = float(arr[3])
                
                
                line = next(f)
                arr = line.split()
                pot_height = float(arr[3])

                while(True):
                    line = next(f)
                    arr = line.split()
                    pot_cv_num = int(arr[3])
                    line = next(f)
                    arr = line.split()
                    pot_cv_coord = float(arr[2])
                    line = next(f)
                    arr = line.split()
                    pot_width = float(arr[3])

                    gau_pot = GaussianPotentialModel(pot_height, pot_width, pot_cv_num, pot_cv_coord)
                    local_gau_pots.append(gau_pot)
                    
                    try:
                        line = next(f)

                        if ('GAUSSIAN BIAS POTENTIAL' in line or len(line.strip())==0):
                            self.gau_pots.append(local_gau_pots)
                            break
                    except StopIteration:
                        self.gau_pots.append(local_gau_pots)
                        break

    def loadData_fes(self, filename):
        with open(filename, 'r') as f:
            self.fes = []

            coords = []
            values = []
            for line in f:
                if line.strip()[0] == '#':
                    if (len(coords) > 0):
                        fes_step = FESPotential()

                        fes_step.coords = coords
                        fes_step.values = values
                        self.fes.append(fes_step)
                        coords = []
                        values = []
                else:
                    # print(line)
                    arr = line.split()
                    coords.append(list(map(float, arr[0:-1])))
                    values.append(float(arr[-1]))
            if (len(coords) > 0):
                fes_step = FESPotential()

                fes_step.coords = coords
                fes_step.values = values
                self.fes.append(fes_step)
                coords = []
                values = []

        # if len(self.coords) > 0:
        #     self.plot_fes(0)
        #     self.slider.setMaximum(len(self.coords)-1)




class MTDTool(qtw.QMainWindow):
    
    def __init__(self):
        super().__init__()
        
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

        self.cv_coord_tab = CVCoordWindow()
        self.tab.addTab(self.cv_coord_tab, 'CV')

        self.cv_height_tab = CVHeightWindow()
        self.tab.addTab(self.cv_height_tab, 'CV Height')

        self.time_property = TimeCoursePropertyWindow()
        self.tab.addTab(self.time_property, 'Time vs Property')

        self.layout.addWidget(self.tab)

        self.setCentralWidget(main_widget)


        loadFesAct = qtw.QAction(QIcon(), '&Load Data', self)
        loadFesAct.triggered.connect(self.load_data)
        fileMenu.addAction(loadFesAct)
        fileMenu.addAction(exitAct)

        # finalizing
        self.setMinimumSize(800, 600)
        self.setWindowTitle('DCDFTBK MetaDynamics Tool')    
        self.show()

    def load_data(self):
        dialog = qtw.QFileDialog(self, 'Load Data Folder', pathlib.os.curdir)
        folder = dialog.getExistingDirectory()
        
        if folder != "":
            self.folder = pathlib.Path(folder)
            self.model = MetaDynamicsResultModel()
            self.model.parseFolder(self.folder)

            if (len(self.model.gau_pots) == len(self.model.fes) and len(self.model.gau_pots) > 0):
                self.time_slider.setEnabled(True)
                self.time_slider.valueChanged.connect(self.update_plots)

                last_step = (len(self.model.gau_pots)-1)

                self.cur_step_label.setText(str(last_step))
                self.total_step_label.setText(str(last_step))

                self.time_slider.setMaximum(last_step)
                self.time_slider.setValue(last_step)
                self.update_plots(last_step)
                
                
            
    def update_plots(self, step):
        self.cur_step_label.setText(str(step))
        self.cur_time_label.setText('{:8.2f} ps'.format((step+1)*self.model.gaussian_interval_time/1000.0))
        self.fes_tab.plot(step, self.model)
        self.cv_coord_tab.plot(step, self.model)
        self.cv_height_tab.plot(step, self.model)
        
        

if __name__ == '__main__':
    
    app = qtw.QApplication(sys.argv)
    tool = MTDTool()
    sys.exit(app.exec_())

