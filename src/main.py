#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 10:17:42 2024

@author: morita
"""

import numpy as np
from core.file import File
#import CifFile
#from libs.cif2cell.uctools import CellData
from core.elements import numbers
from computation.statistics import histogram,gr,total_gr,ncoeff
import matplotlib.pyplot as plt


path = "./data/amorphous_rmc/sio.cfg"
elements = ['Si','O']
f = File.open(path)
atoms = f.getatoms(0,elements)

#path = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
#f = File.open(path)
#atoms = f.getatoms(0)
symbols = ['Si','O']
dr = 0.05

hist = histogram(atoms.norm_positions,atoms.elements,atoms.volume.vectors,
                 dr,symbols=symbols)

ni = [atoms.numbers[symbols[i]] for i in range(len(symbols))]
r, gr = gr(hist,atoms.volume.vectors,ni,dr)

frac = ni/np.sum(ni)
coeff = ncoeff(symbols,frac)

total_gr = total_gr(gr,coeff)
plt.plot(r,total_gr)
#for i in range(3):    
#    plt.plot(r, gr.T[i])
plt.show()

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
#elements = ['Si','O','Si','O','Si','O','Si','O', 'K', 'H', 'K', 'H']
elements = ['Si','O','Si','O','Si','O','Si','O']
#numbers['SI'] = 1
#numbers['O'] = 2
indexes = list(range(len(elements)))
sorted_elements, sorted_indexes = zip(*sorted(zip(elements, indexes), 
                                              key=lambda x: (numbers.get(str.upper(x[0]), x[1])),
                                              reverse=False))

#print(sorted_elements)
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
