#!/usr/bin/env python3

import pyqtgraph as pg
pg.setConfigOption('background', (255,255,255))
pg.setConfigOption('foreground', (0,0,0))
from scipy.signal import find_peaks
import PyQt5.QtWidgets as qtw
from PyQt5.QtGui import QFont, QColor, QRegExpValidator

from PyQt5.QtCore import QStringListModel, QRegExp


import matplotlib.cm as cmap
import skimage.feature

from .models import *
from .gui_utils import *
from .custom_graphic import *
from .utils import *

class OneDTimeWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()
        
        # self.color_map = ['#1f77b4',  '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',  '#7f7f7f', '#bcbd22', '#17becf']
        self.color_map = ['#4363d8', '#e6194B', '#3cb44b', '#ffe119',  '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabebe', 
            '#469990', '#e6beff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#000000']
        self.model = None
        self.line = None
        
        self.time_scale = 'ps'
        self.time_factor = 1000.0

        self.initUI()

    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        #
        #============================================
        #  plotting
        #============================================
        self.canvas = pg.PlotWidget(axisItems={"left": CustomAxis(orientation="left")})
        self.canvas.setBackground('w')
        self.curves = []
        self.curve_pen = pg.mkPen(color=(0,0,255), width=5, style=QtCore.Qt.SolidLine)
        self.time_lines = []
        self.time_line_pen = pg.mkPen(color=(255, 0, 0), width=3, style=QtCore.Qt.DashLine)

        self.axis_labelStyle = {'font-size': '16pt', 'color': 'black'}
        
        font = QFont()
        font.setPointSize(14)
        self.axis_pen = pg.mkPen(color=(0,0,0), width=5)

        self.canvas.getAxis("bottom").tickFont = font
        self.canvas.getAxis("bottom").setStyle(tickLength=-10, tickTextOffset = 10)
        self.canvas.getAxis('bottom').setPen(self.axis_pen)
        self.canvas.getAxis("bottom").enableAutoSIPrefix(False)
        self.canvas.setLabel('bottom', 'Time', units=self.time_scale, **self.axis_labelStyle)

        self.canvas.getAxis("left").setStyle(tickLength=-10, tickTextOffset = 10)
        self.canvas.getAxis("left").tickFont = font     
        self.canvas.getAxis('left').setPen(self.axis_pen)
        self.canvas.getAxis("left").enableAutoSIPrefix(False)
        self.canvas.getAxis('left').setStyle(showValues=True)

        self.canvas.showAxis('right')
        self.canvas.getAxis('right').setStyle(showValues=False)
        self.canvas.getAxis("right").tickFont = font     
        self.canvas.getAxis('right').setPen(self.axis_pen)
        self.canvas.getAxis("right").enableAutoSIPrefix(False)
        self.canvas.showAxis('top')
        self.canvas.getAxis("top").tickFont = font     
        self.canvas.getAxis('top').setStyle(showValues=False)
        self.canvas.getAxis('top').setPen(self.axis_pen)
        
        #============================================
        self.canvas.getAxis('right').setStyle(showValues=True)
        self.canvas2 = pg.ViewBox()
        self.canvas.scene().addItem(self.canvas2)
        self.canvas.getAxis('right').linkToView(self.canvas2)
        self.canvas2.setXLink(self.canvas)
        self.canvas2.setGeometry(self.canvas.getViewBox().sceneBoundingRect())
        self.updateViews()
        self.canvas.getViewBox().sigResized.connect(self.updateViews)
        #============================================

        self.legend = None
        self.layout.addWidget(self.canvas)
    
    def updateViews(self):
        self.canvas2.setGeometry(self.canvas.getViewBox().sceneBoundingRect())
        self.canvas2.linkedViewChanged(self.canvas.getViewBox(), self.canvas2.XAxis)

    def set_model(self, model):
        if self.model == model:
            pass
        else:
            self.model = model
            self.update_plot()
    

    def update_plot(self):
        for item in self.curves:
            self.canvas.removeItem(item) 
            self.canvas2.removeItem(item) 
        self.curves = []
              

        self.times = self.model.get_times()
        xs = [x/self.time_factor for x in self.times]

        total_line_index = 0

        legend_offset = (-80,30)
        # cur_legend_pos = None
        if self.legend is not None:
            legend_offset = self.legend.offset
            # cur_legend_pos = self.legend.pos()
            try:
                self.legend.scene().removeItem(self.legend)
            except Exception as e:
                print (e)
            

        self.legend = MyLegendItem(size=None, offset=legend_offset)
        self.legend.setParentItem(self.canvas.plotItem)
        # if cur_legend_pos is not None:
            # self.legend.setPos(cur_legend_pos)

        nlines = self.model.get_num_lines_leftaxis()
        if (nlines > 0):
            self.canvas.getAxis('left').setStyle(showValues=True)
            for line_index in range(nlines):
                self.canvas.setLabel('left', self.model.get_value_label_leftaxis(), 
                units=self.model.get_value_unit_leftaxis(), **self.axis_labelStyle)

                pen = pg.mkPen(QColor(self.color_map[total_line_index]), width=5, style=QtCore.Qt.SolidLine)
                total_line_index += 1

                
                curve = self.canvas.plot(xs, self.model.get_values_leftaxis(line_index), pen=pen, antialias=True)
                self.curves.append(curve)
                legend = MyItemSample(curve)
                self.legend.addItem(legend, name=self.model.get_name_leftaxis(line_index) )
        else:
            self.canvas.getAxis('left').setStyle(showValues=False)
                    
        nlines = self.model.get_num_lines_rightaxis()
        if nlines > 0:
            self.canvas.getAxis('right').setStyle(showValues=True)
            for line_index in range(nlines):
                self.canvas.setLabel('right', self.model.get_value_label_rightaxis(), units=self.model.get_value_unit_rightaxis(), **self.axis_labelStyle)

                pen = pg.mkPen(QColor(self.color_map[total_line_index]), width=5, style=QtCore.Qt.SolidLine)
                total_line_index += 1

                curve = pg.PlotDataItem(xs, self.model.get_values_rightaxis(line_index), pen=pen, antialias=True)
                self.canvas2.addItem(curve)
                self.curves.append(curve)
                legend = MyItemSample(curve)
                self.legend.addItem(legend, name=self.model.get_name_rightaxis(line_index) )
        else:
            self.canvas.getAxis('right').setStyle(showValues=False)
        



    def update_time(self, step):          
        if self.model is not None:
            for item in self.time_lines:
                self.canvas.removeItem(item)       
            x_pos = self.times[step]/self.time_factor
            line = self.canvas.addLine(x=x_pos, pen=self.time_line_pen)
            self.time_lines.append(line)

class FESWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.color_map = ['#1f77b4',  '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',  '#7f7f7f', '#bcbd22', '#17becf']
        self.initUI()


    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        self.canvas = pg.PlotWidget(axisItems={"left": CustomAxis(orientation="left")})
        self.canvas.setBackground('w')
        self.curve = None
        self.layout.addWidget(self.canvas)

        self.axis_labelStyle = {'font-size': '16pt', 'color': 'black'}
        self.canvas.setLabel('left', 'Free energy', units='kcal/mol', **self.axis_labelStyle )
        font = QFont()
        font.setPointSize(14)

        self.canvas.getAxis("bottom").tickFont = font
        self.canvas.getAxis("left").tickFont = font
        

        pen = pg.mkPen(color=(0,0,0), width=5)

        self.canvas.getAxis("bottom").setStyle(tickLength=-10, tickTextOffset = 10)
        self.canvas.getAxis("left").setStyle(tickLength=-10, tickTextOffset = 10)
        self.canvas.plotItem.getAxis('left').setPen(pen)
        self.canvas.plotItem.getAxis('bottom').setPen(pen)
        self.canvas.getAxis("left").enableAutoSIPrefix(False)

        self.canvas.showAxis("top")
        self.canvas.showAxis("right")
        self.canvas.getAxis("top").setStyle(tickLength=0, tickTextOffset = 10)
        self.canvas.getAxis("right").setStyle(tickLength=0, tickTextOffset = 10)
        self.canvas.plotItem.getAxis('right').setPen(pen)
        self.canvas.plotItem.getAxis('top').setPen(pen)
        self.canvas.getAxis("right").enableAutoSIPrefix(False)
        self.canvas.getAxis('top').setStyle(showValues=False)
        self.canvas.getAxis('right').setStyle(showValues=False)

        self.lines = []


    def find_minima_barriers(self, coords, values):
        maxima, _ = find_peaks(values)
        minima, _ = find_peaks([-x for x in values])

        return (maxima, minima)

    def plot(self, step, model):
        xs, value_raw = model.get_fes_step(step)      
        values = [x*627.5095 for x in value_raw]

        # if 1d
        if (model.get_fes_dimension() == 1):
            maxima, minima = self.find_minima_barriers(xs, values)
            if self.curve is None:
                pen = pg.mkPen(color=QColor(self.color_map[0]), width=5)
                self.curve = self.canvas.plot(xs, values, pen=pen, antialias=True)
                self.canvas.setLabel('bottom', model.get_cvs()[0][0], **self.axis_labelStyle)
            else:
                self.curve.setData(xs, values)
            
            for item in self.lines:
                self.canvas.removeItem(item)

            index = zip(minima,minima[1:])
            for step, ind in enumerate(index):
                ind_min_1, ind_min_2 = ind
                ind_max = maxima[step]
             
                x_min = xs[ind_min_1]
                x_max = xs[ind_min_2]
                # x_width = x_max-x_min
                y_max = values[ind_max]


                line = self.canvas.plot([x_min, x_max], [y_max, y_max], pen = pg.mkPen(color=(255, 128, 0), width=3, style=QtCore.Qt.DashLine))
                self.lines.append(line)
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
            pass
class TwoDFESWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

        self.minima = {}
        self.maxima = {}
        # self.paths = {}
    
    def dijkstra(self, V, start):
        mask = V.mask
        visit_mask = mask.copy() # mask visited cells
        m = np.ones_like(V) * np.inf
        connectivity = [(i,j) for i in [-1, 0, 1] for j in [-1, 0, 1] if (not (i == j == 0))]
        cc = start # current_cell
        m[cc] = 0
        P = {}  # dictionary of predecessors 
        #while (~visit_mask).sum() > 0:
        for _ in range(V.size):
            neighbors = [tuple(e) for e in np.asarray(cc) - connectivity 
                        if e[0] > 0 and e[1] > 0 and e[0] < V.shape[0] and e[1] < V.shape[1]]
            neighbors = [ e for e in neighbors if not visit_mask[e] ]
            tentative_distance = np.asarray([V[e]-V[cc] for e in neighbors])
            for i,e in enumerate(neighbors):
                d = tentative_distance[i] + m[cc]
                if d < m[e]:
                    m[e] = d
                    P[e] = cc
            visit_mask[cc] = True
            m_mask = np.ma.masked_array(m, visit_mask)
            cc = np.unravel_index(m_mask.argmin(), m.shape)
        return m, P

    def shortestPath(self, start, end, P):
        Path = []
        step = end
        while True:
            Path.append(step)
            if step == start: break
            step = P[step]
        Path.reverse()
        return np.asarray(Path)

    def initUI(self):
        self.layout = qtw.QHBoxLayout(self)
        self.setLayout(self.layout)

        self.canvas = pg.PlotItem(axisItems={"left": CustomAxis(orientation="left")})        
        self.curves = []
        self.lines = []

        self.axis_labelStyle = {'font-size': '16pt', 'color': 'black'}
        
        font = QFont()
        font.setPointSize(14)

        self.canvas.getAxis("bottom").tickFont = font
        self.canvas.getAxis("left").tickFont = font
        

        self.pen = pg.mkPen(color=(0,0,0), width=5)

        self.canvas.getAxis("bottom").setStyle(tickLength=-10, tickTextOffset = 10)
        self.canvas.getAxis("left").setStyle(tickLength=-10, tickTextOffset = 10)
        self.canvas.getAxis('left').setPen(self.pen)
        self.canvas.getAxis('bottom').setPen(self.pen)
        self.canvas.getAxis("left").enableAutoSIPrefix(False)

        self.canvas.showAxis("top")
        self.canvas.showAxis("right")
        self.canvas.getAxis("top").setStyle(tickLength=0, tickTextOffset = 10)
        self.canvas.getAxis("right").setStyle(tickLength=0, tickTextOffset = 10)
        self.canvas.getAxis('right').setPen(self.pen)
        self.canvas.getAxis('top').setPen(self.pen)
        self.canvas.getAxis("right").enableAutoSIPrefix(False)
        self.canvas.getAxis('top').setStyle(showValues=False)
        self.canvas.getAxis('right').setStyle(showValues=False)

        self.colorLegendItem = ColorLegendItem(showHistogram=False)
        self.colorLegendItem.setMinimumHeight(60)
        
        
        self.graphicsLayoutWidget = pg.GraphicsLayoutWidget()
        self.graphicsLayoutWidget.addItem(self.canvas, 0, 0)
        self.graphicsLayoutWidget.addItem(self.colorLegendItem, 0, 1)
        self.layout.addWidget(self.graphicsLayoutWidget)

    def plot(self, step, model):
        data, value_raw = model.get_fes_step(step)      
        x, y = data
        z = value_raw*627.5095

        # if 2d
        if (model.get_fes_dimension() == 2):

            self.canvas.setLabel('left', 'CV2 ({})'.format(model.get_cvs()[1][0]), **self.axis_labelStyle )
            self.canvas.setLabel('bottom', 'CV1 ({})'.format(model.get_cvs()[0][0]), **self.axis_labelStyle )
            self.canvas.getAxis('left').setPen(self.pen)
            self.canvas.getAxis('bottom').setPen(self.pen)

            for item in self.curves:
                self.canvas.removeItem(item)
            self.curves = []
            for item in self.lines:
                self.canvas.removeItem(item)
            self.lines = []

            pos, rgba_colors = zip(*cmapToColormap(cmap.gnuplot2))
            pgColormap =  pg.ColorMap(pos, rgba_colors)
            levels=(z.min(), z.max())
            image = pg.ImageItem(z, autoDownsample=True, levels=levels)
            image.setLookupTable(pgColormap.getLookupTable()) 
            self.canvas.addItem(image)
            self.curve = image
            self.colorLegendItem.setImageItem(image)
            self.colorLegendItem.setLevels(levels)
            
            x_width = (x[-1]-x[0])
            y_width = (y[-1]-y[0])

            
            image.scale(x_width/100.0, y_width/100.0)        
            image.translate(y[0]/y_width*100, x[0]/x_width*100)
            
            minima_pos = skimage.feature.peak_local_max(-z)
            if step in self.minima:
                minima = self.minima[step]
            else:
                minima = []
                for c in minima_pos:
                    minima.append( (x[c[0]], y[c[1]], z[c[0],c[1]] ))

                minima.sort(key=lambda x: x[2])
                self.minima[step] = minima   

            if len(minima_pos) >= 2:    
                font = QFont()
                font.setPointSize(14)
                for minimum in minima:
                    text = pg.TextItem('{:.2f}'.format(minimum[2]), color=(255,255,255), anchor=(-0.1,0.5) )
                    text.setPos(minimum[0], minimum[1])
                    text.setFont(font)
                    self.canvas.addItem(text)
                    self.lines.append(text)
                    item = pg.ScatterPlotItem([minimum[0]], [minimum[1]], size=8, pen=pg.mkPen(color=(0,255,0)) )
                    self.canvas.addItem(item)
                    self.lines.append(item)
                    
                if step in self.maxima:
                    maximum = self.maxima[step]
                else:
                    start = (minima_pos[0][0], minima_pos[0][1])
                    end =  (minima_pos[1][0], minima_pos[1][1])
                    
                    V = np.ma.masked_array(z, z>0)
                    D, P = self.dijkstra(V, start)
                    path = self.shortestPath(start, end, P)

                    path_x = []
                    path_y = []
                    maximum = (0, 0, -1E20)
                    for po in path:
                        path_x.append(y[po[0]])
                        path_y.append(x[po[1]])
                        if z[po[0], po[1]] > maximum[2]:
                            maximum = (x[po[0]], y[po[1]], z[po[0], po[1]])
                    self.maxima[step] = maximum
                
                text = pg.TextItem('{:.2f}'.format(maximum[2]), color=(0,0,0), anchor=(-0.1,0.5) )
                text.setPos(maximum[0], maximum[1])
                text.setFont(font)
                self.canvas.addItem(text)
                self.lines.append(text)
                item = pg.ScatterPlotItem([maximum[0]], [maximum[1]], size=8, pen=pg.mkPen(color=(0,255,0)) )
                self.canvas.addItem(item)
                self.lines.append(item)
                
        else:
            pass

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

        self.textbox.setFocus()

class TimeCoursePropertyWindow(qtw.QWidget):
    def __init__(self, model):
        super().__init__()
        self.initUI()
        
        self.result_model = model
        self.property_selection_model = QStringListModel()
        self.selection_dialog = PropertySelectionWidget(self.property_selection_model)
        
        self.selections = set()
        self.property_selection_changed = False
        self.plot_selection_model = TimeCourseDataModel()
        
        self.plot_axis1 = 'Bond Distance'
        self.plot_axis2 = 'None'
        self.change_axes()
        

    def initUI(self):
        self.layout = qtw.QVBoxLayout(self)
        self.setLayout(self.layout)
        
        self.plotWidget = OneDTimeWindow()
        
        self.hbox = qtw.QHBoxLayout()
        self.layout.addLayout(self.hbox)
        self.layout.addWidget(self.plotWidget)

        self.hbox.setAlignment(Qt.AlignLeft)
        self.button_setting = qtw.QPushButton('Setting')
        self.button_setting.setMaximumWidth(100)
        self.hbox.addWidget(self.button_setting)
        self.button_setting.pressed.connect(self.show_selection)
        

    def update_plot(self):
        times = [self.result_model.gaussian_interval_time*(i+1) for i in range(len(self.result_model.get_gau_pots()))]
        self.plot_selection_model.set_times(times)
        self.plotWidget.set_model(self.plot_selection_model)
        

    def change_axes(self):        
        if self.plot_axis1 == 'Bond Distance':
            self.plot_selection_model.set_leftaxis(self.plot_axis1, 'Angstrom')
        elif self.plot_axis1 == 'Bond Angle':
            self.plot_selection_model.set_leftaxis(self.plot_axis1, 'Degree')
        elif self.plot_axis1 == 'Charge':
            self.plot_selection_model.set_leftaxis(self.plot_axis1, 'e')
        elif self.plot_axis1 == 'None':
            self.plot_selection_model.set_leftaxis(None, 'e')

        if self.plot_axis2 == 'Bond Distance':
            self.plot_selection_model.set_rightaxis(self.plot_axis2, 'Angstrom')
        elif self.plot_axis2 == 'Bond Angle':
            self.plot_selection_model.set_rightaxis(self.plot_axis2, 'Degree')
        elif self.plot_axis2 == 'Charge':
            self.plot_selection_model.set_rightaxis(self.plot_axis2, 'e')
        elif self.plot_axis2 == 'None':
            self.plot_selection_model.set_rightaxis(None, 'e')


    def show_selection(self):
        ret = self.selection_dialog.exec()
        if ret == 0:
            return
        
        selected = self.property_selection_model.stringList()
        self.plot_axis1 = self.selection_dialog.left_y.currentText()
        self.plot_axis2 = self.selection_dialog.right_y.currentText()
        
        if self.selection_dialog.axis_changed:
            self.property_selection_changed = True
            self.change_axes()

        # check 

        new_selection = set()
        for item in selected:
            arr = item.split('-')
            if (len(arr) == 1):
                #charge
                new_selection.add(PlottableProperty('Charge', int(arr[0])))
            elif (len(arr) == 2):
                new_selection.add(PlottableProperty('Bond Distance', int(arr[0]), int(arr[1])))
            elif (len(arr) == 3):
                new_selection.add(PlottableProperty('Bond Angle', int(arr[0]), int(arr[1]), int(arr[2])))
            else:
                raise ValueError('wrong values')

        items_for_deletion = self.selections - new_selection
        for item in items_for_deletion:
            self.plot_selection_model.remove_data(str(item))
            self.selection_changed = True
            self.selections.remove(item)

        
        items_for_insertion = new_selection - self.selections
        for item in items_for_insertion:
            self.property_selection_changed = True

            if item.type == 'Charge':
                act = Actions(self, 'Loading Mulliken...', self.result_model.get_charge)
                act.exec()
                act = Actions(self, 'Loading Trajectory...', self.result_model.get_trajectory)
                act.exec()
                traj = self.result_model.get_trajectory()
                charge = self.result_model.get_charge()
                y_data = [ ts[item.index[0]]  for ts in charge]
                name='C({}{})'.format(self.result_model.symbols[item.index[0]], item.index[0])
                self.plot_selection_model.add_data(item.type, str(item), name, y_data)
                self.selections.add(item)
            elif item.type == 'Bond Distance':
                act = Actions(self, 'Loading Trajectory...', self.result_model.get_trajectory)
                act.exec()
                traj = self.result_model.get_trajectory()
                y_data = [ get_distance(ts[item.index[0]], np.array(ts[item.index[1]]), self.result_model.lattice)  for ts in traj]
                name='R({}{}-{}{})'.format(self.result_model.symbols[item.index[0]], item.index[0], self.result_model.symbols[item.index[1]], item.index[1])
                self.plot_selection_model.add_data(item.type, str(item), name, y_data)
                self.selections.add(item)
            elif item.type == 'Bond Angle':
                act = Actions(self, 'Loading Trajectory...', self.result_model.get_trajectory)
                act.exec()
                traj = self.result_model.get_trajectory()
                y_data = [ get_angle(ts[item.index[0]], np.array(ts[item.index[1]]), np.array(ts[item.index[2]]), self.result_model.lattice)  for ts in traj]
                name='A({}{}-{}{}-{}{})'.format(self.result_model.symbols[item.index[0]], item.index[0], self.result_model.symbols[item.index[1]], item.index[1]
                , self.result_model.symbols[item.index[2]], item.index[2])
                self.plot_selection_model.add_data(item.type, str(item), name, y_data)
                self.selections.add(item)
            
        if self.property_selection_changed:
            self.plotWidget.update_plot()
            self.plotWidget.update_time(self.cur_step)
            self.property_selection_changed = False


    def update_time(self, step):  
        self.cur_step = step
        self.plotWidget.update_time(step)