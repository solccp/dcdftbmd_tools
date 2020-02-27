#!/usr/bin/env python3

import numpy as np
import math
import pymatgen as mg
import pymatgen.util.coord as uc
import collections

def get_distance(cart_coord_1, cart_coord_2, lattice=None):
    if lattice is None:
        vec = cart_coord_1 - cart_coord_2
        return np.linalg.norm(vec)
    else:
        latt = mg.Lattice(lattice)
        fc1 = latt.get_fractional_coords(cart_coord_1)
        fc2 = latt.get_fractional_coords(cart_coord_2)
        res = uc.pbc_shortest_vectors(latt, fc1, fc2, return_d2=True)
        return math.sqrt(res[1][0][0])

def get_angle(cart_coord_1, cart_coord_2, cart_coord_3, lattice=None):
    if (lattice == None):
        vec_1 = np.array(cart_coord_1) - np.array(cart_coord_2)
        vec_2 = np.array(cart_coord_3) - np.array(cart_coord_2)
        return math.degrees(math.acos(np.dot(vec_1, vec_2) / (np.linalg.norm(vec_1) * np.linalg.norm(vec_2))))
    else:
        latt = mg.Lattice(lattice)
        fc1 = latt.get_fractional_coords(cart_coord_1)
        fc2 = latt.get_fractional_coords(cart_coord_2)
        fc3 = latt.get_fractional_coords(cart_coord_3)
        vec_1 = uc.pbc_shortest_vectors(latt, fc2, fc1)[0][0]
        vec_2 = uc.pbc_shortest_vectors(latt, fc2, fc3)[0][0]
        return math.degrees(math.acos(np.dot(vec_1, vec_2) / (np.linalg.norm(vec_1) * np.linalg.norm(vec_2))))

def cmapToColormap(cmap, nTicks=16):
    """
    Converts a Matplotlib cmap to pyqtgraphs colormaps. No dependency on matplotlib.
    Parameters:
    *cmap*: Cmap object. Imported from matplotlib.cm.*
    *nTicks*: Number of ticks to create when dict of functions is used. Otherwise unused.
    """

    # Case #1: a dictionary with 'red'/'green'/'blue' values as list of ranges (e.g. 'jet')
    # The parameter 'cmap' is a 'matplotlib.colors.LinearSegmentedColormap' instance ...
    if hasattr(cmap, '_segmentdata'):
        colordata = getattr(cmap, '_segmentdata')
        if ('red' in colordata) and isinstance(colordata['red'], collections.Sequence):

            # collect the color ranges from all channels into one dict to get unique indices
            posDict = {}
            for idx, channel in enumerate(('red', 'green', 'blue')):
                for colorRange in colordata[channel]:
                    posDict.setdefault(colorRange[0], [-1, -1, -1])[idx] = colorRange[2]

            indexList = list(posDict.keys())
            indexList.sort()
            # interpolate missing values (== -1)
            for channel in range(3):  # R,G,B
                startIdx = indexList[0]
                emptyIdx = []
                for curIdx in indexList:
                    if posDict[curIdx][channel] == -1:
                        emptyIdx.append(curIdx)
                    elif curIdx != indexList[0]:
                        for eIdx in emptyIdx:
                            rPos = (eIdx - startIdx) / (curIdx - startIdx)
                            vStart = posDict[startIdx][channel]
                            vRange = (posDict[curIdx][channel] - posDict[startIdx][channel])
                            posDict[eIdx][channel] = rPos * vRange + vStart
                        startIdx = curIdx
                        del emptyIdx[:]
            for channel in range(3):  # R,G,B
                for curIdx in indexList:
                    posDict[curIdx][channel] *= 255

            rgb_list = [[i, posDict[i]] for i in indexList]

        # Case #2: a dictionary with 'red'/'green'/'blue' values as functions (e.g. 'gnuplot')
        elif ('red' in colordata) and isinstance(colordata['red'], collections.Callable):
            indices = np.linspace(0., 1., nTicks)
            luts = [np.clip(np.array(colordata[rgb](indices), dtype=np.float), 0, 1) * 255 \
                    for rgb in ('red', 'green', 'blue')]
            rgb_list = zip(indices, list(zip(*luts)))

    # If the parameter 'cmap' is a 'matplotlib.colors.ListedColormap' instance, with the attributes 'colors' and 'N'
    elif hasattr(cmap, 'colors') and hasattr(cmap, 'N'):
        colordata = getattr(cmap, 'colors')
        # Case #3: a list with RGB values (e.g. 'seismic')
        if len(colordata[0]) == 3:
            indices = np.linspace(0., 1., len(colordata))
            scaledRgbTuples = [(rgbTuple[0] * 255, rgbTuple[1] * 255, rgbTuple[2] * 255) for rgbTuple in colordata]
            rgb_list = zip(indices, scaledRgbTuples)

        # Case #4: a list of tuples with positions and RGB-values (e.g. 'terrain')
        # -> this section is probably not needed anymore!?
        elif len(colordata[0]) == 2:
            rgb_list = [(idx, (vals[0] * 255, vals[1] * 255, vals[2] * 255)) for idx, vals in colordata]

    # Case #X: unknown format or datatype was the wrong object type
    else:
        raise ValueError("[cmapToColormap] Unknown cmap format or not a cmap!")
    
    # Convert the RGB float values to RGBA integer values
    return list([(pos, (int(r), int(g), int(b), 255)) for pos, (r, g, b) in rgb_list])


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
        return hash(self.__str__())

    def __str__(self):
        return '{}'.format('-'.join(map(str, self.index)))

    def __repr__(self):
        return self.__str__()