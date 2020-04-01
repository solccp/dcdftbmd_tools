#!/usr/bin/env python3

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dcdftbmd_tools.md_gui.models import *
import pathlib

import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np
import pathlib
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 16})
matplotlib.use('QT5Agg')
matplotlib.rcParams["figure.dpi"] = 200.0
matplotlib.rcParams["figure.figsize"] = [10, 6]

def find_minima_barriers(coords, values):
    maxima, _ = find_peaks(values)
    minima, _ = find_peaks([-x for x in values])

    return (maxima, minima)



raw_models = []

paths = sys.argv[1:]
abs_paths = [os.path.abspath(p) for p in paths]

biaspot_ranges = []
biaspot_lens = []


for path in paths:
    model = MetaDynamicsResultModel()
    model.set_base_folder(pathlib.Path(path))
    biaspot_lens.append(len(model.get_gau_pots()))
    min_x, max_x = model.gaussian_coord_range[0]
    width = max_x - min_x
    min_x -= 0.5*width
    max_x += 0.5*width 
    biaspot_ranges.append( (min_x, max_x) )
    raw_models.append(model)



if len(set(biaspot_lens)) > 1:
    print('Warning: Number of bias potentials is consistent: ', (biaspot_lens))
else:
    print('Number of bias potentials: ', biaspot_lens[0])

# print(biaspot_lens)
step = min(biaspot_lens)-1

all_min = min([r[0] for r in biaspot_ranges])
all_max = max([r[1] for r in biaspot_ranges])
print(all_min, all_max)

grids = np.linspace(all_min, all_max, 150)

avg_fes = np.zeros(len(grids), dtype=np.float)
fes = []
for ind, model in enumerate(raw_models):
    step = len(model.get_gau_pots())-1
    _, temp_values = model.get_fes_step(step, (all_min, all_max))
    values = np.array([x*627.5095 for x in temp_values])
    fes.append(values)
    avg_fes += values
    
    plt.plot(grids, values, '--', color='grey')

avg_fes /= float(len(fes))
plt.plot(grids, avg_fes, label='Averaged', color='blue', lw=3)

plt.ylabel('Free energy [kcal/mol]')
plt.xlabel(raw_models[0].get_cvs()[0][0])
plt.legend(frameon=False)

maxima, minima = find_minima_barriers(grids, avg_fes)
index = zip(minima,minima[1:])
for step, ind in enumerate(index):
    ind_min_1, ind_min_2 = ind
    ind_max = maxima[step]

    x_min = grids[ind_min_1]
    x_max = grids[ind_min_2]
    x_width = x_max-x_min
    y_max = avg_fes[ind_max]
    plt.plot([x_min, x_max], [y_max, y_max], '--' )

    min_1_x_pos = grids[ind_min_1]
    min_1_y_pos = avg_fes[ind_min_1]
    min_2_x_pos = grids[ind_min_2]
    min_2_y_pos = avg_fes[ind_min_2]

    energy_1 = y_max - min_1_y_pos
    energy_2 = y_max - min_2_y_pos

    plt.annotate('', xy=(min_1_x_pos, min_1_y_pos), xytext=(min_1_x_pos, y_max), arrowprops=dict(arrowstyle='<|-|>'))

    plt.text(min_1_x_pos+x_width*0.01, (min_1_y_pos+y_max)/2.0, '{:.2f}'.format(energy_1), fontsize=12, 
    horizontalalignment='left', verticalalignment='center', bbox=dict(boxstyle='square,pad=0.1', fc='w', ec='none'))

    plt.annotate('', xy=(min_2_x_pos, min_2_y_pos), xytext=(min_2_x_pos, y_max), arrowprops=dict(arrowstyle='<|-|>'))
    plt.text(min_2_x_pos-x_width*0.01, (min_2_y_pos+y_max)/2.0, '{:.2f}'.format(energy_2), fontsize=12, 
    horizontalalignment='right', verticalalignment='center', bbox=dict(boxstyle='square,pad=0.1', fc='w', ec='none'))           
    
plt.show()



