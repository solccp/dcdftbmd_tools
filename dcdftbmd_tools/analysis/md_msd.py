#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import pymatgen as mg
import pymatgen.util.coord as cu

def run(opt):
    traj_file = opt.traject

    if opt.tube_axis > 0:
        print('Tube axis:', opt.tube_axis, file=sys.stderr)

    index2print = []

    if (opt.elem):
        print("", file=sys.stderr)
        print('Elements to be analyzed: ', ", ".join(opt.elem), file=sys.stderr)
    elif (opt.index):
        index2print = opt.index
    elif (opt.index_file):
        with open(opt.index_file, 'r') as f:
            for line in f:
                line_list = list(map(int, line.split()))
                index2print += line_list
    index2print = sorted(index2print)

    symbols = []
    symbol_list = []
    nat_list = []
    lattice = []

    with open(opt.latt, 'r') as f:
        for line in f:
            if 'TV' in line:
                lattice.append(list(map(float, line.split()[1:4])))

    box = mg.Lattice(lattice)

    #read the first geometry from the trajectory file
    with open(traj_file, 'r') as f:
        nat = int(f.readline())
        info = f.readline()
        arr = info.split()
        t1 = float(arr[3])
        for i in range(0, nat):
            line = f.readline()
            arr = line.split()
            symbols.append(arr[0])
            if (opt.elem and (arr[0] in opt.elem )):
                index2print.append(i)
            if (len(symbol_list) == 0):
                symbol_list.append(arr[0])
                nat_list.append(1)
            else:
                if (arr[0] == symbol_list[-1]):
                    nat_list[-1] += 1
                else:
                    symbol_list.append(arr[0])
                    nat_list.append(1)
        nat = int(f.readline())
        info = f.readline()
        arr = info.split()
        t2 = float(arr[3])
        dt = t2-t1

    index_list = [index2print[x:x+5] for x in range(0, len(index2print),5)]
    index2print_set = set(index2print)

    print ("", file=sys.stderr)
    print ('Indexes to be analyzed:', file=sys.stderr)
    num = 0
    for li in index_list:
        print ('{:3d}:'.format(num*5+1), ('  {:5d}'*len(li)).format(*li), file=sys.stderr)
        num+=1
    print ("", file=sys.stderr)

    if (not opt.index_file):
        with open('index_auto.dat', 'w') as f:
            for li in index_list:
                print (('  {:5d}'*len(li)).format(*li), file=f)

    first = True
    first_coords = []

    step = 1
    # fout = sys.stdout
    fout = open(opt.output, 'w')

    if opt.tube_axis == 0:
        print ('#Time(ps) msd(A^2/ps)', file=fout)
    else:
        print ('#Time(ps) msd msd(parallel) msd(perpendicular)', file=fout)

    # times = []
    # msds = []
    # msds_x = []
    # msds_yz = []

    min_half_box = min(box._lengths)*0.5
    min_half_box = min_half_box*min_half_box

    nat_print = len(index2print)
    coords = np.zeros((nat_print, 3))
    first_coords = np.zeros((nat_print, 3))
    previous_coords = np.zeros((nat_print, 3))

    with open(traj_file, 'r') as f:
        for line in f:

            info = f.readline()
            arr = info.split()
            ind = 0
        
            for i in range(0, nat):
                line = f.readline()
                if step >= opt.start:
                    if (i in index2print_set):
                        arr = line.split()
                        for j in range(3):
                            coords[ind][j] = float(arr[j+1])
                        ind += 1
            
            if step < opt.start:
                step += 1
                continue
            if (first):
                first_coords = coords.copy()
                first = False
                previous_coords = coords.copy()
            else:
                summ = 0.0
                summ_x = 0.0
                summ_yz = 0.0

                for i in range(nat_print):
                    pre_site = previous_coords[i]
                    cur_site = coords[i]
                    
                    dv = cur_site - pre_site
                    
                    dis = dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]
                    
                    if (dis >= min_half_box):
                        frac_pre = box.get_fractional_coords(pre_site)
                        frac_cur = box.get_fractional_coords(cur_site)

                        dv = cu.pbc_shortest_vectors(box, frac_pre, frac_cur)[0][0]
                        new_site = pre_site + dv
                        coords[i] = new_site
                        dv = new_site - first_coords[i]
                    else:
                        dv = cur_site - first_coords[i]
    
                    if opt.tube_axis == 0:
                        lx = dv[0]*dv[0]
                        lyz = (dv[1]*dv[1]+dv[2]*dv[2]) 
                    elif opt.tube_axis == 1:
                        lx = dv[0]*dv[0]
                        lyz = (dv[1]*dv[1]+dv[2]*dv[2]) 
                    elif opt.tube_axis == 2:
                        lx = dv[1]*dv[1]
                        lyz = (dv[0]*dv[0]+dv[2]*dv[2]) 
                    elif opt.tube_axis == 3:
                        lx = dv[2]*dv[2]
                        lyz = (dv[0]*dv[0]+dv[1]*dv[1]) 
                    
                    summ_x += lx
                    summ_yz += lyz 
                    summ += (lx+lyz) 
                
                msd = (summ/nat_print)
                time = dt*step
                
                if opt.tube_axis > 0:
                    msd_x = (summ_x/nat_print)
                    msd_yz = (summ_yz/nat_print)
                    print('{:10.4f} {:20.12f} {:20.12f} {:20.12f}'.format(time/1000.0, msd, msd_x, msd_yz), file=fout)
                else:
                    print('{:10.4f} {:20.12f}'.format(time/1000.0, msd), file=fout)
                
                step += 1
                tmp = previous_coords
                previous_coords = coords
                coords = tmp



