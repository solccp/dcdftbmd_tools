#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib
import matplotlib.pylab as plt
import math


import PyQt5.QtWidgets as qtw
from PyQt5.QtGui import QIcon,  QRegExpValidator, QFont
from PyQt5.QtCore import Qt, QStringListModel, QRegExp, QThread, pyqtSignal
from PyQt5 import QtCore

import pyqtgraph as pg
from scipy.signal import find_peaks

import pathlib
matplotlib.use('QT5Agg')
matplotlib.rcParams["figure.dpi"] = 200.0
matplotlib.rcParams["figure.figsize"] = [10, 6]

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class External(QThread):
    finished = pyqtSignal()

    def __init__(self, callback, parent=None):
        super().__init__(parent=parent)
        self.callback = callback
    
        

    def run(self):
        self.callback()  
        self.finished.emit()

class Actions(qtw.QDialog):
    def __init__(self, callback):
        super().__init__()
        self.initUI()
        self.calc = External(callback)
        self.calc.finished.connect(self.onFinished)
        self.calc.start()
        
    def initUI(self):
        self.setWindowTitle('Progress Bar')
        self.progress = qtw.QProgressBar(self)
        self.progress.setGeometry(0, 0, 300, 25)
        self.progress.setMinimum(0)
        self.progress.setMaximum(0)
        self.progress.setValue(0)
        self.show()

    def onFinished(self):
        self.done(1)
        

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

class CustomAxis(pg.AxisItem):
    @property
    def nudge(self):
        if not hasattr(self, "_nudge"):
            self._nudge = 0
        return self._nudge

    @nudge.setter
    def nudge(self, nudge):
        self._nudge = nudge
        s = self.size()
        # call resizeEvent indirectly
        self.resize(s + QtCore.QSizeF(1, 1))
        self.resize(s)

    def resizeEvent(self, ev=None):
        # s = self.size()

        ## Set the position of the label
        nudge = self.nudge
        br = self.label.boundingRect()
        p = QtCore.QPointF(0, 0)
        if self.orientation == "left":
            p.setY(int(self.size().height() / 2 + br.width() / 2))
            p.setX(-nudge)
        elif self.orientation == "right":
            p.setY(int(self.size().height() / 2 + br.width() / 2))
            p.setX(int(self.size().width() - br.height() + nudge))
        elif self.orientation == "top":
            p.setY(-nudge)
            p.setX(int(self.size().width() / 2.0 - br.width() / 2.0))
        elif self.orientation == "bottom":
            p.setX(int(self.size().width() / 2.0 - br.width() / 2.0))
            p.setY(int(self.size().height() - br.height() + nudge))
        self.label.setPos(p)
        # print(p)
        self.picture = None

    def _updateWidth(self):
        
        if not self.isVisible():
            w = 0
        else:
            if self.fixedWidth is None:
                if not self.style['showValues']:
                    w = 0
                elif self.style['autoExpandTextSpace'] is True:
                    w = self.textWidth
                else:
                    w = self.style['tickTextWidth']
                w += self.style['tickTextOffset'][0] if self.style['showValues'] else 0
                w += max(0, self.style['tickLength'])
                if self.label.isVisible():
                    w += self.label.boundingRect().height() * 1.2 ## Modified
            else:
                w = self.fixedWidth
        
        self.setMaximumWidth(w)
        self.setMinimumWidth(w)
        self.picture = None
        
        


class FESWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        self.canvas = pg.PlotWidget(axisItems={"left": CustomAxis(orientation="left")})
        self.canvas.setBackground('w')
        # self.figure = plt.figure()
        # self.canvas = FigureCanvas(self.figure)
        # self.axis = self.figure.add_subplot(1,1,1)
        # self.figure.subplots_adjust(top=0.950, bottom=0.15, left=0.1, right=0.95)
        # self.toolbar = NavigationToolbar(self.canvas, self)
        # self.layout.addWidget(self.toolbar)
        self.curve = None
        self.layout.addWidget(self.canvas)

        self.axis_labelStyle = {'font-size': '16pt', 'color': 'black'}
        # self.axis_tickStyle = {'font-size': '16pt', 'color': 'black', 'tickLength':10}
        self.canvas.setLabel('left', 'Free energy', units='kcal/mol', **self.axis_labelStyle )
        font = QFont()
        font.setPointSize(14)

        self.canvas.getAxis("bottom").tickFont = font
        self.canvas.getAxis("left").tickFont = font
        self.canvas.getAxis("bottom").setStyle(tickLength=10, tickTextOffset = 10)
        self.canvas.getAxis("left").setStyle(tickLength=10, tickTextOffset = 10)
        # self.canvas.getAxis("bottom").setStyle()
        # self.canvas.getAxis("left").setStyle(tickTextOffset = 20)

        pen = pg.mkPen(color=(0,0,0), width=5)
        self.canvas.plotItem.getAxis('left').setPen(pen)
        self.canvas.plotItem.getAxis('bottom').setPen(pen)
        # self.canvas.getAxis("left").setWidth(120)
        self.canvas.getAxis("left").enableAutoSIPrefix(False)

        self.lines = []
        # self.canvas.re


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
            
            if self.curve is None:
                pen = pg.mkPen(color=(0, 0, 255), width=5)
                self.curve = self.canvas.plot(xs, values, pen=pen, antialias=True)
                self.canvas.setLabel('bottom', model.cvs[0], **self.axis_labelStyle)
            else:
                self.curve.setData(xs, values)
            
            # self.axis.cla()
            # self.axis.plot(xs, values)

            # self.axis.set_xlabel(model.cvs[0])
            # self.axis.set_ylabel('Free energy [kcal/mol]')

            for item in self.lines:
                self.canvas.removeItem(item)

            index = zip(minima,minima[1:])
            for step, ind in enumerate(index):
                ind_min_1, ind_min_2 = ind
                ind_max = maxima[step]
             
                x_min = xs[ind_min_1]
                x_max = xs[ind_min_2]
                x_width = x_max-x_min
                y_max = values[ind_max]


                line = self.canvas.plot([x_min, x_max], [y_max, y_max], pen = pg.mkPen(color=(255, 128, 0), width=3, style=QtCore.Qt.DashLine))
                
                # line = self.canvas.addLine(y=y_max)
                self.lines.append(line)
                # print(line.angle)
                min_1_x_pos = xs[ind_min_1]
                min_1_y_pos = values[ind_min_1]
                min_2_x_pos = xs[ind_min_2]
                min_2_y_pos = values[ind_min_2]

                energy_1 = y_max - min_1_y_pos
                energy_2 = y_max - min_2_y_pos

                line = self.canvas.plot([min_1_x_pos, min_1_x_pos], [min_1_y_pos, y_max], pen = pg.mkPen(color=(0, 0, 0), width=3, style=QtCore.Qt.SolidLine))
                self.lines.append(line)
                
                label1 = pg.TextItem('{:.2f}'.format(energy_1), anchor=(0.0, 0.5), color=(0,0,0))
                label1.setPos(min_1_x_pos, (y_max+min_1_y_pos)*0.5)
                self.lines.append(label1)
                self.canvas.addItem(label1)

                line = self.canvas.plot([min_2_x_pos, min_2_x_pos], [min_2_y_pos, y_max], pen = pg.mkPen(color=(0, 0, 0), width=3, style=QtCore.Qt.SolidLine))
                self.lines.append(line)
                label2 = pg.TextItem('{:.2f}'.format(energy_2), anchor=(1.0, 0.5), color=(0,0,0))
                label2.setPos(min_2_x_pos, (y_max+min_2_y_pos)*0.5)
                self.lines.append(label2)
                self.canvas.addItem(label2)
                               
        else:
            raise RuntimeError('More than 1D FES plotting is not implemented')
class PropertySelectionWidget(qtw.QDialog):
    def __init__(self, model, parent=None, flags=Qt.WindowFlags()):
        super().__init__(parent=parent, flags=flags)
        
        self.model = model
        self.axis_changed = False
        self.initUI()


    def remove_item(self):
        ind = self.tableview.currentIndex()
        # index = self.model.index(ind)        
        self.model.removeRow(ind.row())

    def add_item(self):
        item = self.textbox.text()
        
        pos = 0
        status, value, _ = self.validator.validate(item, pos)
        if status == 2:
            data = self.model.stringList()
            data.append(value)
            self.model.setStringList(data)
            

    def set_axis_changed(self, value):
        self.axis_changed = True

    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        self.hbox = qtw.QGridLayout()
        self.layout.addLayout(self.hbox)
        self.setWindowTitle('Property Selector')

        self.setMinimumSize(600, 480)

        self.hbox.addWidget(qtw.QLabel('Left Y Axis:'), 0, 0)
        self.hbox.addWidget(qtw.QLabel('Right Y Axis:'), 1, 0)

        self.axis_model = QStringListModel(['Bond Distance', 'Bond Angle', 'Charge', 'None'])
        
        self.left_y = qtw.QComboBox()
        self.left_y.setInsertPolicy(qtw.QComboBox.NoInsert)
        self.left_y.setModel(self.axis_model)
        self.hbox.addWidget(self.left_y, 0, 1)

        self.right_y = qtw.QComboBox()
        self.right_y.setInsertPolicy(qtw.QComboBox.NoInsert)
        self.right_y.setModel(self.axis_model)
        self.right_y.setCurrentText('None')
        self.left_y.currentTextChanged.connect(self.set_axis_changed)

        self.right_y.currentTextChanged.connect(self.set_axis_changed)
        self.hbox.addWidget(self.right_y, 1, 1)
        tip = qtw.QLabel('Tip: Single value for charge, 1-2 for bond distance between index 1 and 2, and 1-2-3 for bond angle')
        tip.setWordWrap(True)
        self.hbox.addWidget(tip, 2, 0, 2, 2)

        self.textbox = qtw.QLineEdit()
        rx = QRegExp(r"\d+|\d+[-]\d+|\d+[-]\d+[-]\d+")
        self.validator = QRegExpValidator(rx)

        self.hbox.addWidget(self.textbox, 4,0,1,2)

        self.button_add = qtw.QPushButton('Add')
        self.button_add.pressed.connect(self.add_item)
        self.button_remove = qtw.QPushButton('Remove')
        self.button_remove.pressed.connect(self.remove_item)

        self.hbox.addWidget(self.button_add, 5, 0)
        self.hbox.addWidget(self.button_remove, 5, 1)
        self.button_ok = qtw.QPushButton('OK')
        
        
        self.tableview = qtw.QListView()
        self.tableview.setModel(self.model)
        self.button_ok.pressed.connect(self.accept)
        self.layout.addWidget(self.tableview)
        self.button_ok.setMaximumWidth(100)
        self.layout.addWidget(self.button_ok)

class PlottableProperty:
    def __init__(self, type_, *index):
        super().__init__()
        self.type = type_
        self.index = index
    
    def __eq__(self, other):
        return other and self.type == other.type and self.index == other.index

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash('{}-{}'.format(self.type, '-'.join(map(str, self.index))))

    def __str__(self):
        return '{}-{}'.format(self.type, '-'.join(map(str, self.index)))
    def __repr__(self):
        return self.__str__()

def getAngle(vec_A, vec_B, vec_C):
    vec_1 = np.array(vec_A) - np.array(vec_B)
    vec_2 = np.array(vec_C) - np.array(vec_B)
    return math.acos(np.dot(vec_1, vec_2) / (np.linalg.norm(vec_1) * np.linalg.norm(vec_2)))

class TimeCoursePropertyWindow(qtw.QWidget):
    def __init__(self, model):
        super().__init__()
        self.initUI()
        self.model = model
        self.line = None
        self.traj = []
        self.mull = []
        self.selection_model = QStringListModel()
        self.selection_dialog = PropertySelectionWidget(self.selection_model)
        
        self.selections = {}
        self.selected_plot_objects = {}
        self.selection_changed = False
        self.ylim = [0, 0]
        self.color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.figure.subplots_adjust(top=0.950, bottom=0.15, left=0.1, right=0.95)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.axis = self.figure.add_subplot(1,1,1)
        self.axis2 = self.axis.twinx()
        self.hbox = qtw.QHBoxLayout()
        self.layout.addLayout(self.hbox)
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)

        self.hbox.setAlignment(Qt.AlignLeft)
        self.button_setting = qtw.QPushButton('Setting')
        self.button_setting.setMaximumWidth(100)
        self.hbox.addWidget(self.button_setting)
        self.button_setting.pressed.connect(self.show_selection)
       

    def plot(self, step, model):
        self.cur_step = step
        if (self.line):
            self.axis.lines.pop()
            # self.axis.lines.remove(self.line)
            self.line = None

        if self.selection_changed == False:
            pass
        else:
            # self.model = model

            self.axis.cla()

            axis1_type = 'None'
            if self.plot_axis1 == 'Bond Distance':
                self.axis.set_ylabel(r'Bond Distance [$\AA$]')
                axis1_type = 'R'
            elif self.plot_axis1 == 'Bond Angle':
                self.axis.set_ylabel(r'Bond Angle [degree]')
                axis1_type = 'A'
            elif self.plot_axis1 == 'Charge':
                self.axis.set_ylabel(r'Charge')
                axis1_type = 'C'
            
            axis2_type = 'None'
            if self.plot_axis2 == 'Bond Distance':
                self.axis2.set_ylabel(r'Bond Distance [$\AA$]')
                axis2_type = 'R'
            elif self.plot_axis2 == 'Bond Angle':
                self.axis2.set_ylabel(r'Bond Angle [degree]')
                axis2_type = 'A'
            elif self.plot_axis2 == 'Charge':
                self.axis2.set_ylabel(r'Charge')
                axis2_type = 'C'
            
            ind = 0
            for name, values in self.selections.items():
                name_info = str(name).split('-')
                name_type = name_info[0]
                if name_info[0] == 'bondlength':
                    name_type = 'R'
                elif  name_info[0] == 'bondangle':
                    name_type = 'A'
                elif name_info[0] == 'charge':
                    name_type = 'C'
                
                xs = [self.model.gaussian_interval_time*x/1000.0 for x in range(len(values))]
                if name_type == axis1_type:
                    obj = self.axis.plot(xs, values, label='{}({})'.format(name_type, '-'.join(name_info[1:])), color=self.color_cycle[ind])
                    self.selected_plot_objects[name] = obj
                    ind += 1
                elif name_type == axis2_type:
                    obj = self.axis2.plot(xs, values, label='{}({})'.format(name_type, '-'.join(name_info[1:])), color=self.color_cycle[ind])
                    self.axis.plot(np.nan, label='{}({})'.format(name_type, '-'.join(name_info[1:])), color=self.color_cycle[ind])
                    self.selected_plot_objects[name] = obj
                    ind += 1
            
            self.ylim = self.axis.get_ylim()

            self.axis.set_xlabel('Time [ps]')
            # self.axis.set_ylabel('Gaussian Height of CV1 [kcal/mol]')
            
            
            # labs = [l.get_label() for l in self.selected_plot_objects.values()]
            # self.axis.legend(self.selected_plot_objects, labs)
            self.axis.legend()
            # self.axis2.legend()
            self.selection_changed = False
        self.axis.set_xlim(left=0.0)
        
        
    
        x_pos = self.model.gaussian_interval_time*(step+1)/1000.0
        
        self.line = self.axis.plot([x_pos, x_pos], [self.ylim[0], self.ylim[1]], '--', color='red')
        self.canvas.draw()


    def show_selection(self):
        self.selection_dialog.exec()
        
        selected = self.selection_model.stringList()
        self.plot_axis1 = self.selection_dialog.left_y.currentText()
        self.plot_axis2 = self.selection_dialog.right_y.currentText()
        if self.selection_dialog.axis_changed:
            self.selection_changed = True

        # check 

        new_selection = set()
        for item in selected:
            arr = item.split('-')
            if (len(arr) == 1):
                #charge
                new_selection.add(PlottableProperty('charge', int(arr[0])))
            elif (len(arr) == 2):
                new_selection.add(PlottableProperty('bondlength', int(arr[0]), int(arr[1])))
            elif (len(arr) == 3):
                new_selection.add(PlottableProperty('bondangle', int(arr[0]), int(arr[1]), int(arr[2])))
            else:
                raise ValueError('wrong values')

        items_for_deletion = self.selections.keys() - new_selection
        for item in items_for_deletion:
            self.selected_plot_objects.popitem(item)
            self.selections.popitem(item)
            self.selection_changed = True

        
        items_for_insertion = new_selection - self.selections.keys()
        for item in items_for_insertion:
            self.selection_changed = True
            if item.type == 'charge':
                if len(self.mull) == 0:
                    act = Actions(self.load_charge)
                    act.exec()
                ys = [x[item.index[0]] for x in self.mull]
                self.selections[item] = ys
            elif item.type == 'bondlength':
                if len(self.traj) == 0:
                    act = Actions(self.load_trajectory)
                    act.exec()
                ys = [ np.linalg.norm(np.array(ts[item.index[0]][1])-np.array(ts[item.index[1]][1])) for ts in self.traj]
                self.selections[item] = ys
            elif item.type == 'bondangle':
                if len(self.traj) == 0:
                    act = Actions(self.load_trajectory)
                    act.exec()
                ys = [ getAngle(ts[item.index[0]][1], ts[item.index[1]][1], ts[item.index[2]][1]) for ts in self.traj]
                self.selections[item] = ys
            
        if self.selection_changed :
            self.plot(self.cur_step, self.model)
        # print(self.selections)
                            
    def load_charge(self):
        filename = self.model.folderName.joinpath('dftb.mull')
        if not os.path.exists(filename):
            filename = self.model.folderName.joinpath('mulliken')
        with open(filename, 'r') as f:
            for line in f:
                arr = line.split()
                norb = int(arr[0])
                nat = int(arr[1])
                charges = np.zeros(nat)
                title = next(f)
                arr = title.split()
                cur_step = int(arr[9])
                if cur_step % self.model.gaussian_interval_step == 0 :
                    for i in range(norb):
                        line = next(f)
                        arr = line.split()
                        charges[int(arr[0])-1] += float(arr[3])
                    self.mull.append(charges.tolist())
                else:
                    for i in range(norb):
                        line = next(f)

        
    def load_trajectory(self):
        filename = self.model.folderName.joinpath('dftb.traj')
        if not os.path.exists(filename):
            filename = self.model.folderName.joinpath('traject')
        with open(filename, 'r') as f:
            for line in f:
                coords = []
                nat = int(line)
                title = next(f)
                arr = title.split()
                cur_step = int(arr[9])
                if cur_step % self.model.gaussian_interval_step == 0 :
                    for i in range(nat):
                        line = next(f)
                        arr = line.split()
                        coords.append([arr[0], list(map(float, arr[1:4]))])
                    self.traj.append(coords)
                else:
                    for i in range(nat):
                        line = next(f)


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
        self.folderName = folderName
        self.loadData_metacvs(folderName.joinpath('metacv.dat'))
        self.loadData_fes(folderName.joinpath('fes.dat'))
        self.loadData_biaspot(folderName.joinpath('biaspot'))

    def loadData_metacvs(self, fileName):
        with open(fileName, 'r') as f:
            line = next(f)
            if ('RATIONALCOORDINATIONNUMBER' in line):
                self.cvs.append('Coordination Numbers')
            elif ('BONDDISTANCE') in line:
                self.cvs.append('Bond Distance')



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

        self.cv_coord_tab = CVCoordWindow()
        self.tab.addTab(self.cv_coord_tab, 'CV')

        self.cv_height_tab = CVHeightWindow()
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
        self.show()

    def load_data(self, folder):
        self.folder = pathlib.Path(folder)
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
    def clean_data(self):
        self.model = MetaDynamicsResultModel()
        self.time_slider.setEnabled(False)

    def load_data_dialog(self):
        dialog = qtw.QFileDialog(self, 'Load Data Folder', pathlib.os.curdir)
        folder = dialog.getExistingDirectory()
        
        if folder != "":
            self.load_data(folder)
                
                
            
    def update_plots(self, step):
        self.cur_step_label.setText(str(step))
        self.cur_time_label.setText('{:8.2f} ps'.format((step+1)*self.model.gaussian_interval_time/1000.0))
        self.fes_tab.plot(step, self.model)
        self.cv_coord_tab.plot(step, self.model)
        self.cv_height_tab.plot(step, self.model)
        self.time_property.plot(step, self.model)
        
        

if __name__ == '__main__':
    
    app = qtw.QApplication(sys.argv)
    tool = MTDTool()
    try:
        tool.load_data(os.curdir)
    except Exception as e:
        print(e)
        tool.clean_data()
    sys.exit(app.exec_())

