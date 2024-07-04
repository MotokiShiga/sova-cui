# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from sova.core.file import File
from sova.computation.structure_factor import neighbor
import matplotlib.pyplot as plt
import numpy as np

path = "../data/amorphous_rmc/sio.cfg"
f = File.open(path)
elements = ['Si','O']
atoms = f.getatoms(0, elements)

si = [i for i, s in enumerate(atoms.elements) if s == 'Si']
ox = [i for i, s in enumerate(atoms.elements) if s == 'O']

neis = []
hist = np.zeros(10, dtype='int')
for i in si:
    nei, dis = neighbor(atoms, i, 2.8)
    num = 0
    for n in nei:
        if n in ox:
            num += 1
    hist[num] += 1

num = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
hist = np.array(hist)
plt.bar(num, hist)
plt.xlabel('neighbors')
plt.ylabel('number of count')
plt.show()
