# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from sova.core.file import File
from sova.computation.structure_factor import polyhedra
import matplotlib.pyplot as plt
import numpy as np

path = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()

_polyhedra = polyhedra(atoms,center='Si',around='O',rmax=2.0)

values = []
for poly in _polyhedra:
    if poly.q is not None:
        values.append(poly.q)

y, x = np.histogram(values, bins=20, range=[0.9,1.1])
plt.bar((x[:-1]+x[1:])*0.5, y, width=0.8*(x[1]-x[0]))
plt.xlim(0.9,1.1)
plt.ylim(0, None)
plt.show()

