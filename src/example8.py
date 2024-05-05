# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from core.file import File
from computation.statistics import triplets
import matplotlib.pyplot as plt

path = "./data/crystal/si222.cif"
f = File.open(path)
atoms = f.getatoms()

rmax = [2.8]
angle, hist = triplets(atoms,r=rmax)

# show graph
w = angle[1]-angle[0]
fig = plt.figure(figsize=(12, 6)) 
fig.suptitle("{} : TRIPLETS - bond angle correlations".format(path))
ax = fig.add_subplot()
ax.bar(angle, hist[0], width=w*0.8, label=atoms.trios[0])
ax.set_xlim(0.0,180.0)
#ax.set_ylim(0,800)
ax.legend()
plt.show()

# write output
print('Output from TRIPLETS - bond angle correlations')
print()
print('Maximum r values : ')
s = ''
for r in rmax:
    s += "%17.8F" % r
print(s)
print('\n')

for i, trio in enumerate(atoms.trios):
    print('{0}\n'.format(trio))
    print('PLOTS')
    for x, y in zip(angle, hist[i]): 
        print("%17.8F" % x + "%17.8F" % y)
    print()

