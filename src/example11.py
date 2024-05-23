# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from core.file import File
from computation.statistics import rings
import matplotlib.pyplot as plt
import numpy as np

path = "./data/crystal/sio2_alpha_cristobalite222.cif"
f = File.open(path)
atoms = f.getatoms()

ring = rings(atoms)

#print(len(ring.rings))
#print(ring.rings[0])

#index = [0, 4, 25, 31, 24, 28, 1, 7]
#for i in index:
#    print(atoms.elements[i], atoms.positions[i])

"""
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
"""