#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 10:17:42 2024

@author: morita
"""

from core.file import File
import CifFile
from libs.cif2cell.uctools import CellData
from core.elements import numbers
from computation.pdf import histogram

"""
path = "./data/amorphous_rmc/sio.cfg"
elements = ['Si','O']
f = File.open(path)
atoms = f.getatoms(0,elements)
#print(atoms.grid)
"""

path = "./data/amorphous_md/a_Si_speed1e11K.xyz"
f = File.open(path)
atoms = f.getatoms(0)

print(atoms.norm_positions, atoms.elements)

"""
#path = "./data/pyMolDyn/structure_c.xyz"

path = "./data/a_Si.xyz"
f = File.open(path)
atoms = f.getatoms(0)

path = 'data/pyMolDyn/hexagonal.xyz'
f = File.open(path)
atoms = f.getatoms(0)

path = "./data/AlPO-17.cif"
f = File.open(path)
atoms = f.getatoms(0)

cf = CifFile.ReadCif(path, grammar='1.1')
cb = cf.get(cf.keys()[0])

# Get cell data
cell_data = CellData()
cell_data.getFromCIF(cb)
#print(cell_data.a, cell_data.c)
#print(cell_data.gamma)

"""

"""

elements = ['Si','O','Si','O','Si','O','Si','O', 'K', 'H', 'K', 'H']
indexes = list(range(len(elements)))

sorted_elements, sorted_indexes = zip(*sorted(zip(elements, indexes), 
                                              key=lambda x: (numbers.get(str.upper(x[0]), x[1])),
                                              reverse=True))

#print(sorted_data)

#print(elements)
#print(indexes)
#sorted_elements, sorted_indexes = zip(*sorted(zip(elements, indexes)), key=lambda x: numbers.get(str.upper(x)))

#sorted_elements = sorted(elements, xx, key=lambda x: numbers.get(str.upper(x)))

print(sorted_elements)
print(sorted_indexes)

import copy

atomic_numbers = copy.copy(numbers)

print(id(atomic_numbers))
print(id(numbers))

elements = ['K','Al','Si','O']

sorted_elements = sorted(elements, key=lambda x: atomic_numbers.get(str.upper(x)))

print(sorted_elements)
"""
