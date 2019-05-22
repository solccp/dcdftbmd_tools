#!/usr/bin/env python3

import re

class OutputParser:
    def __init__(self, options={}):
        self.program = 'unknown'
        self.version_arch = 'unknown'
        self.version_number = 'unknown'
        self.version_date = 'unknown'
        self.begin_time = 'unknown'
        self.end_time = 'unknown'
        self.keywords = ""
        self.options = options
        pass

    def _parse_geom_opt(self, enumerator):
        for line_no, line in enumerator:
            if 'Molecular coordinate [Angstrom]' in line:
                next(enumerator)
                next(enumerator)
                next(enumerator)
                next(enumerator)
                syms = []
                coords = []
                for i in range(self.no_atoms):
                    _, line = next(enumerator)
                    arr = line.split()
                    sym = arr[0]
                    coord = list(map(float, arr[1:4]))
                    syms.append(sym)
                    coords.append(coord)
                self.last_geom = (syms, coords)
            if 'Lattice vector [Angstrom]' in line:
                next(enumerator)
                next(enumerator)
                next(enumerator)
                next(enumerator)
                syms = []
                coords = []
                for i in range(3):
                    _, line = next(enumerator)
                    arr = line.split()
                    sym = arr[0]
                    coord = list(map(float, arr[1:4]))
                    syms.append(sym)
                    coords.append(coord)
                self.last_lattice = (syms, coords)

    def _parse_md(self, enumerator):
        pass

    def parse(self, filename='dftb.out'):
        with open(filename, 'r') as fin:
            enumerator = enumerate(fin)
            for line_no, line in enumerator:
                if 'DDDDDD      CCCCCC  DDDDDD    FFFFFFFF TTTTTTTTT BBBBBBB   MM     MM DDDDD' in line and line_no < 10:
                    self.program = 'dcdftbmd'
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    line_no, line = next(enumerator)                 
                    m = re.search('\s*(?P<arch>.*) VERSION (?P<ver_num>[0-9.]+) \((?P<ver_date>.*)\)$', line)
                    if m is not None:
                        self.version_arch = m['arch']
                        self.version_number = m['ver_num']
                        self.version_date = m['ver_date']
                elif '*        ######   #####         ######  ####### ####### ######          #   #        *' in line and line_no < 10:
                    self.program = 'dcdftbk'
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    next(enumerator)
                    line_no, line = next(enumerator)             
                    m = re.search('\s*\*\s*(?P<arch>.*) Version (?P<ver_num>[0-9.]+) \((?P<ver_date>.*)\)', line)
                    if m is not None:
                        self.version_arch = m['arch']
                        self.version_number = m['ver_num']
                        self.version_date = m['ver_date']
                m = re.match('\s+Execution of .* begun (?P<begin_date>.*)', line)
                if m:
                    self.begin_time = m['begin_date']
                    next(enumerator)
                    next(enumerator)
                    titles = []
                    for line_no, line in enumerator:
                        if list(set(line.strip())) == ['-']:
                            break
                        titles.append(line.strip())
                    self.title = '\n'.join(titles)
                
                m = re.match('  (?P<section>[A-Z]+.*[A-Za-z]+)\s+=\s+(?P<value>True|False)', line)
                if m:
                    keywords = [line.rstrip()]
                    for line_no, line in enumerator:
                        if 'Total number of basis set shells' in line:
                            self.no_basis_functions_shell = int(line.split()[7])
                            break
                        keywords.append(line.rstrip())
                    self.keywords = '\n'.join(keywords)
                if 'Number of basis' in line:
                    self.no_basis_functions = int(line.split()[5])
                if 'Number of electrons' in line:
                    self.no_electrons = int(line.split()[4])
                if 'Charge of system' in line:
                    self.charge = int(line.split()[4])
                if 'Spin multiplicity' in line:
                    self.spin_multiplicity = int(line.split()[3])
                if 'Number of occupied orbitals' in line:
                    self.no_occ_orbitals = int(line.split()[5])
                if 'Total number of atoms' in line:
                    self.no_atoms = int(line.split()[5])
                m = re.match('\s*\*\*\* Start (?P<run_type>.*) \*\*\*\s*', line)
                if m:
                    self.run_type = m['run_type']
                    if self.run_type == 'geometry optimization':
                        self._parse_geom_opt(enumerator)
                    elif self.run_type == 'molecular dynamics':
                        self._parse_md(enumerator)
                    else:
                        print('Parsing ', self.run_type , ' is not implemented')
                    return 


    