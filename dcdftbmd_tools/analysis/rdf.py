#!/usr/bin/env python3
import sys
import matplotlib
matplotlib.use('PS')
matplotlib.rc('text', usetex=True)
#matplotlib.use('Qt5agg')
import matplotlib.pyplot as plt
import MDAnalysis
import MDAnalysis.analysis.rdf
import numpy as np
import pymatgen as mg
import argparse
import importlib
import builtins
import math


def rdf():
    sys.path.append('.')

    parser = argparse.ArgumentParser(description='RDF calulater')
    parser.add_argument('-s', '--start', default=1, type=int)
    parser.add_argument('-e', '--endstep', default=-1, type=int)
    parser.add_argument('-r', '--endrange', default=10.0, type=float)

    # parser.add_argument('-e', '--end', default=-1, type=int)
    parser.add_argument('-l', '--latt', type=str, required=True, help='Lattice')

    parser.add_argument('filename', metavar='filename', type=str)
    parser.add_argument('-o', '--output', type=str)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--input', type=str, default='atomgroup', help='AtomGroup definition')
    group.add_argument('-p', '--pair', type=str, help='Atom pair "O,O" for O-O')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-n', '--nbins', default=75, type=int, help='Number of bins')
    group.add_argument('-d', '--resolution', default=0.01, type=float, help='Resolution for bins')

    opt = parser.parse_args(sys.argv[1:])



    outfile = sys.stdout
    make_plot = False
    if (opt.output):
        outfile = open('{}.out'.format(opt.output), 'w')
        make_plot = True

    lattice = []

    with open(opt.latt, 'r') as f:
        for line in f:
            if 'TV' in line:
                lattice.append(list(map(float, line.split()[1:4])))

    box = mg.Lattice(lattice)
    (lx, ly, lz), (a, b, c) = box.lengths_and_angles


    u = MDAnalysis.Universe(opt.filename, format='XYZ')
    u.dimensions = np.array([lx, ly, lz, a, b, c], dtype=np.float32)

    builtins.u = u

    ag1 = None
    ag2 = None
    if opt.pair is not None:
        elems = opt.pair.split(',')
        if (len(elems))>1:
            ag1 = u.select_atoms(f'name {elems[0]}')
            ag2 = u.select_atoms(f'name {elems[1]}')
            with open('atomgroups.py', 'w') as f:
                print(f"ag1 = u.select_atoms('name {elems[0]}')", file=f)
                print(f"ag2 = u.select_atoms('name {elems[1]}')", file=f)
    else:
        atomgroups = importlib.import_module(opt.input)

        ag1 = atomgroups.ag1
        ag2 = atomgroups.ag2
    # ag2 = u.select_atoms('name Pt')

    bins = opt.nbins
    if (opt.resolution):
        nbins = int((opt.endrange - 0.05)/opt.resolution)

    rdf = MDAnalysis.analysis.rdf.InterRDF(ag1, ag2, nbins=nbins, range=(0.05, opt.endrange), verbose=True)
    rdf.run(start=opt.start, stop=opt.endstep, verbose=True)

    # rdf = MDAnalysis.analysis.rdf.InterRDF_s(u, [[ag1,ag2]], nbins=nbins, range=(0.05, opt.endrange))
    # rdf.run(start=opt.start, stop=opt.endstep, verbose=True)

    vol = np.power(rdf.edges[1:], 3) - np.power(rdf.edges[:-1], 3)
    vol *= 4/3.0 * np.pi   
    box_vol = box.volume
    pair_den = box_vol / (len(ag1.atoms) * len(ag2.atoms))
    integral = rdf.rdf / (pair_den / vol)


    #transform 
    #fact = 0.5
    #final_rdf = np.zeros(len(rdf.rdf), dtype=np.double)
    #for i, r1 in enumerate(rdf.bins):
    #    for r, v in zip(rdf.bins, rdf.rdf):
    #        final_rdf[i] += 1.0 + (v-1.0)*math.exp(-(abs(r-r1)-1))**2


    total = 0.0
    step = 0
    for r, v, d in zip(rdf.bins, rdf.rdf, integral):
        step += 1
        total += d/len(ag1.atoms)
        if (step == 1):
            continue
        print('{:8.4f} {:8.4f} {:8.4f}'.format(r, v, total), file=outfile)

    outfile.close()
    if (make_plot):
        plt.plot(rdf.bins, rdf.rdf, 'k-')
        plt.xlabel('r[$\\AA$]')    
        plt.ylabel('g(r)')
        plt.savefig('{}.eps'.format(opt.output), dpi=1200, format='eps')
