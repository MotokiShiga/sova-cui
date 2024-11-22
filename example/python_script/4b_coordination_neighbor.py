from sovapy.core.file import File
from sovapy.computation.structure_factor import neighbor
import matplotlib.pyplot as plt
import numpy as np

import sovapy
print('sovapy ver: ', sovapy.__version__)

# Load structural information from a xyz file
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K.xyz"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

# Setting of maximum bond lengths for atomic element pairs
# Set -1 to pairs for which you don't want to build bonds.
bond_lengths = {('Si', 'O') : 2.0, ('Si', 'Si') : -1, ('O', 'O') : -1}

# Indices of Si and O atoms
ids_si = [i for i, s in enumerate(atoms.elements) if s == 'Si']
ids_o = [i for i, s in enumerate(atoms.elements) if s == 'O']

# Initial value of maximum coordination number
# cn_max is increased if the maximamum is larger
cn_max = 1 
# Calculate the histogram of coordination number around Si atom
hist_cn_si = np.zeros(cn_max+1, dtype='int')
for i in ids_si:
    #This computation does NOT change chemical bonds.
    nei, dis = neighbor(atoms, i, bond_lengths)
    num = 0
    for n in nei:
        if n in ids_o:
            num += 1
    if num>=cn_max:
        hist_cn_si = np.r_[hist_cn_si,np.zeros(num-cn_max+1, dtype='int')]
        cn_max = num+1
    hist_cn_si[num] += 1

# Plot distribution of coordination numbers
cood_num = np.arange(cn_max+1)
plt.bar(cood_num, hist_cn_si)
plt.xlabel('Coordination number of Si atom')
plt.ylabel('Counts')
plt.xticks(cood_num)
plt.show()


# Initial value of maximum coordination number
# cn_max is increased if the maximamum is larger
cn_max = 1 
# Calculate the histogram of coordination number around O atom
hist_cn_o = np.zeros(cn_max+1, dtype='int')
for i in ids_o:
    nei, dis = neighbor(atoms, i, bond_lengths)
    num = 0
    for n in nei:
        if n in ids_si:
            num += 1
    if num>=cn_max:
        hist_cn_o = np.r_[hist_cn_o,np.zeros(num-cn_max+1, dtype='int')]
        cn_max = num+1
    hist_cn_o[num] += 1

# Plot distribution of coordination numbers
cood_num = np.arange(cn_max+1)
plt.bar(cood_num, hist_cn_o)
plt.xlabel('Coordination number of O atom')
plt.ylabel('Counts')
plt.xticks(cood_num)
plt.show()