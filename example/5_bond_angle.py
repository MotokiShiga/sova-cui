# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita and M. Shiga
"""

from sova.core.file import File
from sova.computation.structure_factor import triplets
import matplotlib.pyplot as plt

# Load structural information from a xyz file
structure_file = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

print("Atom symbols:", atoms.symbols)

# Thresholds of distance between atom pairs to enumerate neighbors
# Order of pairs are UNKNOWN!
rmax = [1.32, 1.77, 2.22]
# Calculate bond angle distributions
angle, hist = triplets(atoms,r=rmax)

#Plot bond angle distributions
w = angle[1]-angle[0]
fig = plt.figure(figsize=(12, 5)) 
fig.suptitle("{} : TRIPLETS - bond angle correlations".format(structure_file))
for i, trio in enumerate(atoms.trios):
    ax = fig.add_subplot(2, 3, i+1)
    ax.bar(angle, hist[i], width=w*0.8, label=trio)
    ax.set_xlim(0.0,180.0)
    ax.set_ylim(0,1.0)
    ax.set_xlabel('Angle (Degree)')
    ax.set_ylabel('Probability')
    ax.legend()
plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(hspace=0.3)
plt.show()

# # Standard output 
# print('\nOutput from TRIPLETS - bond angle correlations')
# print()
# print('Maximum r values : ')
# s = ''
# for r in rmax:
#     s += "%17.8F" % r
# print(s)
# print('\n')
# for i, trio in enumerate(atoms.trios):
#     print('{0}\n'.format(trio))
#     print('Angle (degree), Probability ')
#     for x, y in zip(angle, hist[i]): 
#         print("%17.8F" % x + "%17.8F" % y)
#     print()
