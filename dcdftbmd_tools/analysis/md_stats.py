#!/usr/bin/env python3

import sys
import argparse
import statistics
import numpy as np


def md_stats():
    parser = argparse.ArgumentParser(description='Stats for MD properties')
    parser.add_argument('-s', '--start', default=0, type=int)
    parser.add_argument('-e', '--end', default=-1, type=int)
    parser.add_argument('-p', '--property', choices=['temp', 'mdenergy', 'kenergy', 'potenergy'], default='temp')
    parser.add_argument('-m', '--merged', action='store_true')
    parser.add_argument('-f', '--figure', action='store_true')
    parser.add_argument('filename', metavar='filename', type=str)
    parser.set_defaults(merged=False)
    opt = parser.parse_args(sys.argv[1:])

    if (not opt.filename):
        parser.print_usage()
        sys.exit(1)

    data = []
    times = []
    if (opt.merged):
        with open(opt.filename, 'r') as f:
            for line in f:
                if (line.startswith('#')):
                    continue
                if opt.property == 'temp':
                    value = float(line.split()[1])
                    data.append(value)
                elif opt.property == 'mdenergy':
                    value = float(line.split()[4])
                    data.append(value)
                elif opt.property == 'kenergy':
                    value = float(line.split()[3])
                    data.append(value)
                elif opt.property == 'potenergy':
                    value = float(line.split()[2])
                    data.append(value)
    else:
        with open(opt.filename, 'r') as f:
            for line in f:
                if "THIS RUN'S STEP NO.=" in line:
                    time = float(line.split()[3])
                    times.append(time)
                if opt.property == 'temp' and 'TEMPERATURE' in line:
                    value = float(line.split()[2])
                    data.append(value)
                elif opt.property == 'mdenergy' and 'TOTAL MD ENERGY' in line:
                    value = float(line.split()[4])
                    data.append(value)
                elif opt.property == 'kenergy' and 'KINETIC ENERGY' in line:
                    value = float(line.split()[3])
                    data.append(value)
                elif opt.property == 'potenergy' and ('Final' in line and 'DFTB' in line and 'Energy' in line) :
                    value = float(line.split()[4])
                    data.append(value)
                #    data = list(map(float, [line.split()[0] for line in f]))

    real_data = data[opt.start:opt.end]
    real_times = np.array(times[opt.start:opt.end])/1000.0
    times_ps = np.array(times)/1000.0
    pmax = max(real_data)
    pmin = min(real_data)

    print ()
    stddev = statistics.stdev(real_data)
    if 'energy' in opt.property:
        print ('{:<6s} {:15s} {:15s} {:15s} {:15s}'.format('N', 'Mean', 'Min', 'Max', 'StdDev(Kcal/mol)'))
        stddev *= 627.5095
    else:
        print ('{:<6s} {:15s} {:15s} {:15s} {:15s}'.format('N', 'Mean', 'Min', 'Max', 'StdDev(K)'))

    mean = statistics.mean(real_data)
    print ('-'*70)
    print('{:<6d} {:<15.8f} {:<15.8f} {:<15.8f} {:<15.8f}'.format( len(real_data), mean, pmin, pmax, stddev))

    print ()
    #print(statistics.median(real_data), statistics.pstdev(real_data))
    if opt.figure:
        import pyqtgraph as pg
        pw = pg.PlotWindow('MD Statistics', background='w')

        font = pg.QtGui.QFont()
        font.setPointSize(14)
        axis_pen = pg.mkPen(color=(0,0,0), width=3)

        pw.getAxis("bottom").tickFont = font
        pw.getAxis("bottom").setStyle(tickLength=-10, tickTextOffset = 10)
        pw.getAxis('bottom').setPen(axis_pen)
        pw.getAxis("bottom").enableAutoSIPrefix(False)
        axis_labelStyle = {'font-size': '16pt', 'color': 'black'}
        pw.setLabel('bottom', 'Time', units='ps', **axis_labelStyle)

        pw.getAxis("left").setStyle(tickLength=-10, tickTextOffset = 10)
        pw.getAxis("left").tickFont = font     
        pw.getAxis('left').setPen(axis_pen)
        pw.getAxis("left").enableAutoSIPrefix(False)
        pw.getAxis('left').setStyle(showValues=True)
       
        # pw.addLine(x=times_ps[opt.start], pen=pg.mkPen(color=(0,255,0), width=2))
        final_stddev = stddev
        if opt.property == 'temp':
            pw.getAxis('left').setWidth(80)
            pw.setLabel('left', 'Temperature', units='K', **axis_labelStyle)
        elif 'energy' in opt.property:
            pw.getAxis('left').setWidth(120)
            pw.setLabel('left', 'Energy', units='Hartree', **axis_labelStyle)
            final_stddev = stddev/627.5095
        
        avg_line = pw.addLine(y=mean, pen=pg.mkPen(color=(255,0,0), width=3))
        upper_std = pw.addLine(y=mean+final_stddev, pen=pg.mkPen(color=(255,128,0), width=2))
        lower_std = pw.addLine(y=mean-final_stddev, pen=pg.mkPen(color=(255,128,0), width=2))
        upper_std.setZValue(-5)
        lower_std.setZValue(-5)
        avg_line.setZValue(0)

        avg_linetext = pg.InfLineLabel(avg_line, '{:.3f}'.format(mean), position=0.1, fill=(255,255,255), anchors = [(0.5, 0.5), (0.5, 1)])
        avg_linetext.setFont(font)
        avg_linetext.setColor('k')
        avg_linetext.setZValue(-3)
        def update_lines(mean, stddev):
            avg_line.setValue(mean)
            lower_std.setValue(mean-stddev)
            upper_std.setValue(mean+stddev)
            avg_linetext.setText('{:.3f}'.format(mean))
        
        pen = pg.mkPen(color=(0,0,255), width=3)
        curve = pw.plot(times_ps, data, pen=pen, antialias=True)
        curve.setZValue(-10)

        lr = pg.LinearRegionItem([times_ps[opt.start],times_ps[-1]], bounds=[0,times_ps[-1]])
        lr.setZValue(-100)
        pw.addItem(lr)
        
        def updateSelection(sel):
            ran = sel.getRegion()
            sel.setRegion((ran[0], times_ps[-1]) )
            import bisect
            start = bisect.bisect(times_ps, ran[0])-1
            real_data = data[start:]
            mean = statistics.mean(real_data)
            stddev = statistics.stdev(real_data)
            update_lines(mean, stddev)

        lr.sigRegionChanged.connect(updateSelection)
        if sys.flags.interactive != 1 or not hasattr(pg.QtCore, 'PYQT_VERSION'):
            pg.QtGui.QApplication.exec_()

        
if __name__ == '__main__':
    md_stats()