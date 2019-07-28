#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pymatgen as mg

def run():

    parser = argparse.ArgumentParser(description='Translate Coordinates')
    parser.add_argument('traject', metavar='traject', type=str)
    parser.add_argument('-l', '--latt', type=str, required=False, help='Lattice')
    parser.add_argument('-s', '--shift', type=str, required=True, help='Shift')

    opt = parser.parse_args()

    pbc = False
    TVs = []
    if opt.latt :
        with open(opt.latt, 'r') as f:
            for line in f:
                if 'TV' in line:
                    TVs.append(list(map(float, line.split()[1:4])))

        if len(TVs) == 3:
            pbc = True
    lattice = mg.Lattice(TVs)
    shifts = list(map(float, opt.shift.split(',')))
    if len(shifts) != 3:
        raise RuntimeError('Shifts should be "Sx,Sy,Sz"')

    with open(opt.traject, 'r') as f:
        while(True):
            line = next(f)
            if (line.strip()==''):
                break
            nat = int(line)
            print(nat)
            title = next(f)
            print(title.strip())
            for _ in range(nat):
                line = next(f)
        
                arr = line.split()
                coords = list(map(float, arr[1:4]))
                if pbc == False:
                    new_coords = [x+y for x, y in zip(coords, shifts)]
                    print('{:<4s} {:20.12f} {:20.12f} {:20.12f}'.format(arr[0], *new_coords))
                else:
                    new_frac_coordinates = [x+y for x, y in zip(lattice.get_fractional_coords(coords), shifts)]
                    final_frac_coordinates = []
                    for x in new_frac_coordinates:
                        if x >= 1.0:
                            final_frac_coordinates.append(x-1.0)
                        elif x<0.0:
                            final_frac_coordinates.append(x+1.0)
                        else:
                            final_frac_coordinates.append(x)
                    new_coords = lattice.get_cartesian_coords(final_frac_coordinates)
                    print('{:<4s} {:20.12f} {:20.12f} {:20.12f}'.format(arr[0], *new_coords))
