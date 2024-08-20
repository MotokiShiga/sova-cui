# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita and M. Shiga
"""

from sova.core.file import File
from sova.computation.structure_factor import neighbor
import matplotlib.pyplot as plt
import numpy as np


# Load structural information from a xyz file
structure_file = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

# Indices of Si and O atoms
ids_si = [i for i, s in enumerate(atoms.elements) if s == 'Si']
ids_o = [i for i, s in enumerate(atoms.elements) if s == 'O']

# Initial value of maximum coordination number
# cn_max is increased if the maximamum is larger
cn_max = 1 
# Calculate the histogram of coordination number around Si atom
hist_cn_si = np.zeros(cn_max+1, dtype='int')
for i in ids_si:
    nei, dis = neighbor(atoms, i, 2.0)
    num = 0
    for n in nei:
        if n in ids_o:
            num += 1
    if num>=cn_max:
        hist_cn_si = np.r_[hist_cn_si,np.zeros(num-cn_max+1, dtype='int')]
        cn_max = num+1
    hist_cn_si[num] += 1

# Plot distribution of corrdination numbers
cood_num = np.arange(cn_max+1)
plt.bar(cood_num, hist_cn_si)
plt.xlabel('The number of neighbors around Si')
plt.ylabel('Counts')
plt.xticks(cood_num)
plt.show()


# Initial value of maximum coordination number
# cn_max is increased if the maximamum is larger
cn_max = 1 
# Calculate the histogram of coordination number around O atom
hist_cn_o = np.zeros(cn_max+1, dtype='int')
for i in ids_o:
    nei, dis = neighbor(atoms, i, 2.0)
    num = 0
    for n in nei:
        if n in ids_si:
            num += 1
    if num>=cn_max:
        hist_cn_o = np.r_[hist_cn_o,np.zeros(num-cn_max+1, dtype='int')]
        cn_max = num+1
    hist_cn_o[num] += 1

# Plot distribution of corrdination numbers
cood_num = np.arange(cn_max+1)
plt.bar(cood_num, hist_cn_o)
plt.xlabel('The number of neighbors around O')
plt.ylabel('Counts')
plt.xticks(cood_num)
plt.show()