#!/usr/bin/env python3

import sys
import argparse
import statistics


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
    real_times = times[opt.start:opt.end]
    pmax = max(real_data)
    pmin = min(real_data)

    print ()
    stddev = statistics.stdev(real_data)
    if 'energy' in opt.property:
        print ('{:<6s} {:15s} {:15s} {:15s} {:15s}'.format('N', 'Mean', 'Min', 'Max', 'StdDev(Kcal/mol)'))
        stddev *= 627.5095
    else:
        print ('{:<6s} {:15s} {:15s} {:15s} {:15s}'.format('N', 'Mean', 'Min', 'Max', 'StdDev(K)'))

    print ('-'*70)
    print('{:<6d} {:<15.8f} {:<15.8f} {:<15.8f} {:<15.8f}'.format( len(real_data), statistics.mean(real_data), pmin, pmax, stddev))

    print ()
    #print(statistics.median(real_data), statistics.pstdev(real_data))
    if opt.figure:
        import matplotlib.pyplot as plt
        plt.plot(real_times, real_data)
        plt.show()
        
