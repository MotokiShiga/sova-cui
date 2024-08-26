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

# Setting of maximum bond lengths for atomic element pairs
# Set -1 to pairs for which you don't want to build bonds.
bond_lengths = {('Si', 'O') : 2.0, ('Si', 'Si') : -1, ('O', 'O') : -1}
atoms.set_bond_lengths(bond_lengths)

# Indices of Si and O atoms
ids_si = np.array([i for i, s in enumerate(atoms.elements) if s == 'Si'])
ids_o = np.array([i for i, s in enumerate(atoms.elements) if s == 'O'])

# list of corrdination numbers
cnums = np.array([len(b) for b in atoms.bonds])

# lists of coordination numbers around Si and O atoms
cnums_si = cnums[ids_si]
cnums_o = cnums[ids_o]

# maxima of coordination numbers
cn_max_si = cnums_si.max()
cn_max_o  = cnums_o.max()

# Calculate the histogram of coordination number around Si atom
hist_cn_si = np.zeros(cn_max_si +1, dtype='int')
for cn in range(cn_max_si+1):
    hist_cn_si[cn] = np.sum(cnums_si==cn)

# Plot distribution of corrdination numbers
cood_num = np.arange(cn_max_si+1)
plt.bar(cood_num, hist_cn_si)
plt.xlabel('The number of neighbors around Si')
plt.ylabel('Counts')
plt.xticks(cood_num)
plt.show()

# Calculate the histogram of coordination number around Si atom
hist_cn_o = np.zeros(cn_max_o +1, dtype='int')
for cn in range(cn_max_o+1):
    hist_cn_o[cn] = np.sum(cnums_o==cn)

# Plot distribution of corrdination numbers
cood_num = np.arange(cn_max_o+1)
plt.bar(cood_num, hist_cn_o)
plt.xlabel('The number of neighbors around O')
plt.ylabel('Counts')
plt.xticks(cood_num)
plt.show()