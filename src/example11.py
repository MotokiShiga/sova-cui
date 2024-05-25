# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from core.file import File
from computation.statistics import rings
from core.gridding import Grid
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
import igraph as ig
from core.molecule import RING

path = "./data/crystal/sio2_beta_cristobalite222.cif"
f = File.open(path)
atoms = f.getatoms()

ring = RING()
ring.set_atoms(atoms)

# GUTTMAN, KING, PRIMITIVE
ring.calculate(ring_type=RING.RingType.GUTTMAN, 
               pair_atom_symbols=[['Si', 'O']], 
               pair_dist_max=[2.0])


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

sx = []; sy = []; sz = []
ox = []; oy = []; oz = []
index = 100
x = []; y = []; z = []
for i in ring.rings[index]:
    _x, _y, _z = atoms.positions[i]
    x.append(_x)
    y.append(_y)
    z.append(_z)
    if atoms.elements[i] == 'Si':
        sx.append(_x)
        sy.append(_y)
        sz.append(_z)
    else:
        ox.append(_x)
        oy.append(_y)
        oz.append(_z)
i = ring.rings[index][0]
_x, _y, _z = atoms.positions[i]
x.append(_x)
y.append(_y)
z.append(_z)

ax.scatter(sx, sy, sz, s=200, c="red", label="Si")
ax.scatter(ox, oy, oz, s=140, c="blue", label="O")
       
line= art3d.Line3D(x, y, z, color='c')
ax.add_line(line)

ax.legend(loc=2)

plt.show()
