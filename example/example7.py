# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from sova.core.file import File
from sova.computation.structure_factor import triplets
import matplotlib.pyplot as plt

path = "../data/sio.cfg"
f = File.open(path)
elements = ['Si','O']
atoms = f.getatoms(0,elements)
#atoms.symbol_order(['Si','O'])
rmax = [1.32, 1.77, 2.22]
angle, hist = triplets(atoms,r=rmax)

w = angle[1]-angle[0]
fig = plt.figure(figsize=(12, 5)) 
fig.suptitle("{} : TRIPLETS - bond angle correlations".format(path))
for i, trio in enumerate(atoms.trios):
    ax = fig.add_subplot(2, 3, i+1)
    ax.bar(angle, hist[i], width=w*0.8, label=trio)
    ax.set_xlim(0.0,180.0)
    ax.set_ylim(0,1.0)
    ax.set_xlabel('angle')
    ax.set_ylabel('distribution')
    ax.legend()
plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(hspace=0.3)
plt.show()

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

"""
path = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()

#nei, dis = neighbor(atoms,1,2.0)
#neighbors(atoms,rmin=[0.,0.,0.], rmax=[4.,2.,3.])

path = "./data/crystal/si222.cif"
#path = "./data/a_Si.xyz"
f = File.open(path)
atoms = f.getatoms()

nei, dis = neighbor(atoms,0,2.8)
print(nei)
#neighbors(atoms,rmin=[0.,0.,0.], rmax=[4.,2.,3.])
"""