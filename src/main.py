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
from computation.statistics import (histogram,gr,total_gr,SQ,total_SQ,total_FQ,
                                    ncoeff,xcoeff,Gr,Tr,Nr)
import matplotlib.pyplot as plt
from core.gridding import volume

#path = "./data/amorphous_rmc/sio.cfg"
#elements = ['Si','O']
#f = File.open(path)
#atoms = f.getatoms(0,elements)

path = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()
dr = 0.05

# calculate histogram
#_, hist = histogram(atoms,dr,symbols=['Si','O'])
_, hist = histogram(atoms,dr)

# calculate g(r)
r, gr = gr(atoms,hist,dr)

# calculate Total g(r)
coeff = ncoeff(atoms.symbols,atoms.frac)
total_gr = total_gr(gr,coeff)

# calculate S(Q)
dq = 0.05
qmin = 0.3
qmax = 25.0
q, sq = SQ(atoms,gr,qmin,qmax,dr,dq)

# calculate Total S(Q)
total_sq = total_SQ(sq,coeff)

# calculate F(Q)
coeff = xcoeff(atoms.symbols,atoms.frac,q)
total_fq = total_FQ(sq,coeff)

# calculate Gr
rho = atoms.rho
_Gr = Gr(r,total_gr,rho)

# calculate Tr
_Tr = Tr(r,total_gr,rho)

# calculate Nr
_Nr = Nr(r,_Tr)

#plt.plot(r,_Nr)
for i in range(3):    
    plt.plot(q, sq.T[i], label=atoms.pairs[i])
plt.legend()
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
