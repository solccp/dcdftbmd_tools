#!/usr/bin/env python3

import numpy as np
import os
import bisect
from copy import deepcopy

class TimeCourseDataModel():
    def __init__(self):
        
        self.times = []

        self.value_label_leftaxis = None
        self.value_label_rightaxis = None

        self.value_unit_leftaxis = None
        self.value_unit_rightaxis = None

        self.data = []

    def set_times(self, times):
        self.times = times

    def get_times(self):
        return self.times

    def set_leftaxis(self, label, unit):
        if label == 'None':
            self.value_label_leftaxis = None
            self.value_unit_leftaxis = None
        else:
            self.value_label_leftaxis = label
            self.value_unit_leftaxis = unit

    def set_rightaxis(self, label, unit):
        if label == 'None':
            self.value_label_rightaxis = None
            self.value_unit_rightaxis = None
        else:
            self.value_label_rightaxis = label
            self.value_unit_rightaxis = unit


    def add_data(self, data_type, id, label, data):
        self.data.append((data_type, id, label, data))

    def remove_data(self, id):
        data = [(t, id_, l, d) for t, id_, l, d in self.data if id_ != id]
        self.data = data

    def get_num_lines_rightaxis(self):
        lines = 0
        for data_type, _, _, _ in self.data:
            if self.value_label_rightaxis == data_type:
                lines += 1
        return lines

    def get_num_lines_leftaxis(self):
        lines = 0
        for data_type, _, _, _ in self.data:
            if self.value_label_leftaxis == data_type:
                lines += 1
        return lines
    
    def get_values_leftaxis(self, index):
        lines = 0
        for data_type, _, _, data_values_y in self.data:
            if self.value_label_leftaxis == data_type:
                if (lines == index):
                    return data_values_y
                lines += 1
        raise ValueError('Out of Index: ', index)
    
    def get_value_label_leftaxis(self):
        return self.value_label_leftaxis

    def get_value_unit_leftaxis(self):
        return self.value_unit_leftaxis
    
    def get_name_leftaxis(self, index):
        lines = 0
        for data_type, _, data_label, _ in self.data:
            if self.value_label_leftaxis == data_type:
                if (lines == index):
                    return data_label
                lines += 1
        raise ValueError('Out of Index: ', index)

    
    def get_values_rightaxis(self, index):
        lines = 0
        for data_type, _, _, data_values_y in self.data:
            if self.value_label_rightaxis == data_type:
                if (lines == index):
                    return data_values_y
                lines += 1
        raise ValueError('Out of Index: ', index)
    
    def get_value_label_rightaxis(self):
        return self.value_label_rightaxis

    def get_value_unit_rightaxis(self):
        return self.value_unit_rightaxis
    
    def get_name_rightaxis(self, index):
        lines = 0
        for data_type, _, data_label, _ in self.data:
            if self.value_label_rightaxis == data_type:
                if (lines == index):
                    return data_label
                lines += 1
        raise ValueError('Out of Index: ', index)

    

class CVHeightAdpater:
    H2kcalmol = 627.5095
    def __init__(self, model):
        self.model = model
        self.times = None 
        self.values = None

    def get_times(self):
        if self.times is None:
            self.times = [self.model.gaussian_interval_time*(i+1) for i in range(len(self.model.get_gau_pots()))]
        return self.times
        
    def get_values_leftaxis(self, index):
        if self.values is None:
            self.values = [x.height*self.H2kcalmol for x in self.model.get_gau_pots()]
        return self.values
    
    def get_value_label_leftaxis(self):
        return 'Gaussian Height'

    def get_value_unit_leftaxis(self):
        return 'kcal/mol'
    
    def get_name_leftaxis(self, index):
        return 'CV'

    def get_num_lines_leftaxis(self):
        return 1

    def get_num_lines_rightaxis(self):
        return 0

class CVCoordAdpater:
    def __init__(self, model):
        self.model = model
        self.times = None 
        self.left_values = None
        self.right_values = None
        
        self.value_unit_leftaxis = self.model.get_cvs()[0][0]
        
        if self.get_num_lines_rightaxis() > 0:
            self.value_unit_rightaxis = self.model.get_cvs()[1][0]    

    def get_times(self):
        if self.times is None:
            self.times = [self.model.gaussian_interval_time*(i+1) for i in range(len(self.model.get_gau_pots()))]
        return self.times
        
    def get_value_label_leftaxis(self):
        return 'CV1'

    def get_value_label_rightaxis(self):
        return 'CV2'

    def get_value_unit_leftaxis(self):
        return self.value_unit_leftaxis

    def get_value_unit_rightaxis(self):
        return self.value_unit_rightaxis
    
    
    def get_name_leftaxis(self, index):
        return 'CV1'
    def get_name_rightaxis(self, index):
        return 'CV2'

    def get_values_leftaxis(self, index):
        if self.left_values is None:
            self.left_values = [x.cv_coords[0] for x in self.model.get_gau_pots()]
        return self.left_values
    
    def get_values_rightaxis(self, index):
        if self.right_values is None:
            self.right_values = [x.cv_coords[1] for x in self.model.get_gau_pots()]
        return self.right_values
    
    def get_num_lines_leftaxis(self):
        return 1

    def get_num_lines_rightaxis(self):
        if self.model.get_fes_dimension() >= 2:
            return 1
        else:
            return 0

class GaussianPotentialModel:
    def __init__(self, height, widths, cv_coords):
        super().__init__()            
        self.height = height
        self.widths = widths
        self.cv_coords = cv_coords
        if len(self.widths) != len(self.cv_coords):
            raise RuntimeError('Dimension is not consistent')
    
    def dimension(self):
        return len(self.widths)

    def value(self, xs):
        if len(xs) != self.dimension():
            raise RuntimeError('Dimension is not consistent')
        exp_value = 0.0
        
        for x, sigma, center in zip(xs, self.widths, self.cv_coords):
            exp_value += (np.power(x - center, 2.) / (2 * np.power(sigma, 2.)))
        value = -self.height*np.exp(-exp_value)
        return value

class MetaDynamicsResultModel:
    
    def __init__(self):
        super().__init__()
        self.gau_pots = None
        self.fes = None
        self.cvs = None
        self.trajectory = None
        self.charge = None
        self.lattice = None
        self.symbols = None
        self.metadata = None
        self.gaussian_interval_step = None
        self.gaussian_interval_time = None
        
        self.default_names = {
            'main_input': ['dftb.inp'],
            'main_output': ['dftb.out'],
            'trajectory': ['traject', 'dftb.traj'],
            'mulliken': ['mulliken', 'dftb.mull'],
            'bias_potential': ['biaspot'],
            'metacv': ['metacv.dat'],
            'fes': ['fes.dat']
        }
    
    def get_fes_dimension(self):
        return self.get_gau_pots()[0].dimension()

    def set_base_folder(self, folderName):
        self.folderName = folderName
        self.gau_pots = None
        self.fes = None
        self.cvs = None
        self.trajectory = None
        self.net_atomic_charge = None
        self.metadata = None
        self.gaussian_interval_step = None
        self.gaussian_interval_time = None
        
    def load_metadata(self, filename):
        with open(filename, 'r') as f:
            for line in f:
                if 'MD         =     True' in line:
                    line = next(f)
                    section = {}
                    for line in f:
                        if len(line.strip()) == 0:
                            break
                        arr = line.split('=')
                        key = arr[0].strip().replace(' ', '_')
                        value = arr[1].strip()
                        section[key] = value
                    if len(section) > 0:
                        if self.metadata is None:
                            self.metadata = {}
                        self.metadata['md'] = section
                elif 'Total number of atoms' in line:
                    self.metadata['Total_number_of_atoms'] = int(line.split()[-1])
                    break


    def get_metaData(self):
        if self.metadata is None:
            for name in self.default_names['main_output']:
                filename = self.folderName.joinpath(name)
                if os.path.exists(filename):
                    self.load_metadata(filename)
                    break
        if self.metadata is None:
            raise RuntimeError('Cannot load file: {}'.format(', '.join(self.default_names['main_output'])))
        return self.metadata
    

    def loadData_metacvs(self, fileName):
        with open(fileName, 'r') as f:
            self.cvs = []
            for line in f:
                if ('RATIONALCOORDINATIONNUMBER' in line):
                    self.cvs.append(('Coordination Numbers', 'coordination_number'))
                elif ('BONDDISTANCE') in line:
                    self.cvs.append(('Bond Distance', 'bond_distance'))
                elif 'ATOMPOINTPLANEDISTANCE' in line:
                    self.cvs.append(('Plane Distance', 'plane_distance'))

    def get_cvs(self):
        if self.cvs is None:
            for name in self.default_names['metacv']:
                filename = self.folderName.joinpath(name)
                if os.path.exists(filename):
                    self.loadData_metacvs(filename)
                    break
        if self.cvs is None:
            raise RuntimeError('Cannot load file: {}'.format(', '.join(self.default_names['metacv'])))
        return self.cvs


    def loadData_biaspot(self, fileName):
        """
        GAUSSIAN BIAS POTENTIAL:        1
        *** AT T=         80.00 FSEC, THIS RUN'S STEP NO.=      80
        Gaussian height     =        0.0010000000 a.u.
        Collective variable =                   1
            Coordinate        =        0.0208250915
            Gaussian width    =        0.1000000000
        """
        min_coord = {}
        max_coord = {}
        first = True
        second = True
        first_step = 0
        first_time = 0.0
        with open(fileName, 'r') as f:
            line = next(f)
            gau_pots = []
            while('GAUSSIAN BIAS POTENTIAL:' in line):
                
                arr = line.split()
                pot_step = int(arr[3])
                title_line = next(f)
                if 'RECOVERED FROM RESTART FILE' not in title_line:
                    if first:
                        first = False
                        arr = title_line.split()
                        first_step = int(arr[9])
                        first_time = float(arr[3])
                    elif second:
                        second = False
                        arr = title_line.split()
                        interval_step = int(arr[9])
                        self.gaussian_interval_step = interval_step - first_step
                        self.gaussian_interval_time = float(arr[3]) - first_time
                        
                
                
                line = next(f)
                arr = line.split()
                pot_height = float(arr[3])      

                pot_widths = []
                pot_coords = []
                line = next(f)
                while(True):
                    arr = line.split()
                    pot_cv_num = int(arr[3])
                    line = next(f)
                    arr = line.split()
                    pot_cv_coord = float(arr[2])

                    if pot_cv_coord > max_coord.get(pot_cv_num-1, -1E20):
                        max_coord[pot_cv_num-1] = pot_cv_coord
                    if pot_cv_coord < min_coord.get(pot_cv_num-1, 1E20):
                        min_coord[pot_cv_num-1] = pot_cv_coord

                    pot_coords.append(pot_cv_coord)

                    line = next(f)
                    arr = line.split()
                    pot_width = float(arr[3])

                    pot_widths.append(pot_width)
                    
                    try:
                        line = next(f)

                        if ('GAUSSIAN BIAS POTENTIAL' in line or len(line.strip())==0):
                            raise StopIteration
                    except StopIteration:
                        pot = GaussianPotentialModel(pot_height, pot_widths, pot_coords)
                        gau_pots.append(pot)
                        break
        if len(gau_pots) > 0:
            self.gau_pots = gau_pots
            ncv = len(min_coord)
            self.gaussian_coord_range = [(min_coord[i], max_coord[i]) for i in range(ncv)]
            if self.gaussian_interval_step is None:
                metadata = self.get_metaData()
                md_timestep = float(metadata['md']['Timestep'])*1E15
                meta_interval = int(metadata['md']['Addition_frequency'])
                self.gaussian_interval_step = meta_interval
                self.gaussian_interval_time = meta_interval*md_timestep

    def get_gau_pots(self):
        names = self.default_names['bias_potential']
        if self.gau_pots is None:
            for name in names:
                filename = self.folderName.joinpath(name)
                if os.path.exists(filename):
                    self.loadData_biaspot(filename)
                    break
        if self.gau_pots is None:
            raise RuntimeError('Cannot load file: {}'.format(', '.join(names)))
        return self.gau_pots

    def get_trajectory(self):
        names = self.default_names['trajectory']
        if self.trajectory is None:
            for name in names:
                filename = self.folderName.joinpath(name)
                if os.path.exists(filename):
                    self.load_trajectory(filename)
                    break
        if self.trajectory is None:
            raise RuntimeError('Cannot load file: {}'.format(', '.join(names)))
        return self.trajectory

    def load_trajectory(self, filename):
        with open(filename, 'r') as f:
            traj = []
            symbols = []
            first_str = True
            for line in f:
                nat = int(line)
                title = next(f)
                arr = title.split()
                cur_step = int(arr[9])
                if cur_step % self.gaussian_interval_step == 0 :
                    coords = np.zeros( (nat, 3) )
                    for i in range(nat):
                        line = next(f)
                        arr = line.split()
                        if first_str:
                            symbols.append(arr[0])
                        coords[i] = list(map(float, arr[1:4]))
                    if cur_step > 0:
                        traj.append(coords)
                else:
                    for i in range(nat):
                        line = next(f)
                if first_str:
                    first_str = False
        if len(traj) > 0:
            self.trajectory = traj
            self.symbols = symbols

        names = self.default_names['main_input']
        lattice = []
        for name in names:
            filename = self.folderName.joinpath(name)
            if os.path.exists:
                with open(filename, 'r') as f:
                    for line in f:
                        if 'TV' in line:
                            arr = line.split()
                            lattice.append(list(map(float, arr[1:4])))
                break
        if (len(lattice) == 3):
            self.lattice = lattice

    def get_charge(self):
        names = self.default_names['mulliken']
        if self.charge is None:
            for name in names:
                filename = self.folderName.joinpath(name)
                if os.path.exists(filename):
                    self.load_charge(filename)
                    break
        if self.charge is None:
            raise RuntimeError('Cannot load file: {}'.format(', '.join(names)))
        return self.charge

    def load_charge(self, filename):
        with open(filename, 'r') as f:
            all_charges = []
            for line in f:
                arr = line.split()
                norb = int(arr[0])
                nat = int(arr[1])
                charges = np.zeros(nat)
                title = next(f)
                arr = title.split()
                cur_step = int(arr[9])
                if cur_step % self.gaussian_interval_step == 0 :
                    for i in range(norb):
                        line = next(f)
                        arr = line.split()
                        charges[int(arr[0])-1] += float(arr[3])
                    if cur_step > 0:
                        # all_charges.append(charges.tolist())
                        all_charges.append(charges)
                else:
                    for i in range(norb):
                        line = next(f)
        if len(all_charges) > 0:
            self.charge = all_charges
    
    def get_fes_step(self, step):

        if self.fes is None:
            n_dimension = self.get_fes_dimension()
            if  n_dimension == 2:
                self.fes = []
                gau_pots = self.get_gau_pots()[:step+1]               
                axis_bounds = []
                
                for cv_index in range(n_dimension):
                    min_x, max_x = self.gaussian_coord_range[cv_index]
                    width = max_x - min_x
                    min_x -= 0.5*width
                    max_x += 0.5*width 
                    axis_bounds.append( (min_x, max_x ))

                npts = 100
                grids_x = np.linspace(axis_bounds[0][0], axis_bounds[0][1], num=npts)
                grids_y = np.linspace(axis_bounds[1][0], axis_bounds[1][1], num=npts)
                
                zz = np.zeros( (len(grids_x), len(grids_y)) )
                
                for pot in gau_pots:
                    min_x_index = bisect.bisect_right(grids_x, pot.cv_coords[0]-pot.widths[0]*3.0)
                    min_y_index = bisect.bisect_right(grids_y, pot.cv_coords[1]-pot.widths[1]*3.0)
                    max_x_index = bisect.bisect_left(grids_x, pot.cv_coords[0]+pot.widths[0]*3.0)
                    max_y_index = bisect.bisect_left(grids_y, pot.cv_coords[1]+pot.widths[1]*3.0)
                    for i in range(min_x_index, max_x_index):
                        for j in range(min_y_index, max_y_index):
                            zz[i,j] += pot.value([grids_x[i],grids_y[j]])
                                
                    self.fes.append( ((grids_x, grids_y),  deepcopy(zz)))
                
            elif n_dimension == 1:
                self.fes = []
                gau_pots = self.get_gau_pots()[:step+1]
                npoints = 150
                cv_index = 0
                min_x, max_x = self.gaussian_coord_range[cv_index]
                width = max_x - min_x
                min_x -= 0.5*width
                max_x += 0.5*width
                xs = np.linspace(min_x, max_x, npoints)
                ys = np.zeros(npoints)

                for pot in gau_pots:
                    for i in range(npoints):
                        ys[i] += pot.value([xs[i]])
                                
                    self.fes.append((xs, np.copy(ys)))
            else:
                raise ValueError('Not implemented')


        return self.fes[step]