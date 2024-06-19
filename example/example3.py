# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from sova.core.file import File
from sova.computation.structure_factor import histogram
import matplotlib.pyplot as plt

# xyz file (no periodic cell)
path = "../data/a_Si.xyz"
f = File.open(path)
atoms = f.getatoms()

print("Periodic : ", atoms.volume.periodic)

dr = 0.05

# calculate histogram
r, hist = histogram(atoms,dr)

fig = plt.figure() 
ax = fig.add_subplot()
ax.bar(r, hist.T[0], width=dr*0.8, label=atoms.pairs[0])
ax.set_xlim(0.0,5.0)
ax.set_ylim(0,600)
ax.legend()
plt.show()
