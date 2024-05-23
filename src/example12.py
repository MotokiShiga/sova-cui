# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from core.file import File
from computation.statistics import rings
from core.gridding import Grid
import matplotlib.pyplot as plt
import numpy as np
import igraph as ig
from core.molecule import RING

#path = "./sio2_alpha_cristobalite222.cif"
#path = "./data/crystal/sio2_alpha_cristobalite222.cif"
path = "./data/crystal/sio2_alpha_cristobalite.cif"
f = File.open(path)
atoms = f.getatoms()

#print(atoms.number)



ring = RING()
ring.set_atoms(atoms)

ring.make_bond_ghost(pair_atom_symbols=[['Si', 'O']], 
                     pair_dist_max=[2.0], 
                     flag_periodicity=True, 
                     p_pair=0.001)
"""
ring.make_bond_periodic(pair_atom_symbols=[['Si', 'O']], 
                     pair_dist_max=[2.0], 
                     flag_periodicity=True, 
                     p_pair=0.001)
"""
#GUTTMAN
#KING
#PRIMITIVE
ring.enumerate_ring(ring_type=RING.RingType.PRIMITIVE, 
                    cutoff_size=24, num_parallel=0) 

import sys
sys.exit()


_positions = []
ghost = []
for i, x in enumerate(atoms.norm_positions):
    _positions.append(x)
    ghost.append(i)
    
for ix in [-1,0,1]:
    for iy in [-1,0,1]:
        for iz in [-1,0,1]:
            if ix == 0 and iy == 0 and iz == 0:
                continue
            for i, x in enumerate(atoms.norm_positions):
                pos = x + np.array([ix*2,iy*2,iz*2])
                _positions.append(pos)
                ghost.append(i)

positions = np.array(_positions)
positions = positions + np.array([3.,3.,3.])
positions /= 3.
positions -= 1.
vectors = np.array(atoms.volume.vectors)

m = np.matrix([[3,0,0],
               [0,3,0],
               [0,0,3]])
vectors *= m
grid = Grid(positions, vectors)

unit_atoms_index = [i for i in range(atoms.number)]

n = len(positions)
dist_max = 2.0
index_atoms = []
index_cells = []
for ic in range(n):        
    grid.neighbours(ic, dist_max, 0, n-1)
    for inei, d, shift in zip(grid.inei, grid.d, np.array(grid.shifts).T):
        if ic < inei:
            index_atoms.append([ic, inei])
            index_cells.append([shift[0],shift[1],shift[2]])

chemical_bond_index = index_atoms
chemical_bond_index = np.array(chemical_bond_index)
atoms_all = [i for i in range(n)]

G = ig.Graph()
G.add_vertices(len(atoms_all))
G.add_edges(chemical_bond_index)

#unit_bond_index

set_rings = set()

bond_index = chemical_bond_index.tolist()
index_cells = np.array(index_cells)

for i in range(chemical_bond_index.shape[0]):
    n0 = chemical_bond_index[i,0]
    n1 = chemical_bond_index[i,1]
    
    if not (n0 in unit_atoms_index) and not (n1 in unit_atoms_index):
        continue
    
    ei = G.get_eid(n0, n1)
    G.delete_edges(ei)
    
    paths = G.get_all_shortest_paths(n0, n1)
    
    for path in paths:
        """
        sum_shift = np.zeros(3,dtype='int8')
        lpath = list(path)
        for ns, ne in zip(lpath, lpath[1:]+lpath[:1]):
            sign = 1
            if ns > ne:
                ns, ne = ne, ns
                sign = -1
            index = bond_index.index([ns,ne])
            shift = index_cells[index]
            sum_shift += sign*shift
        
        # change unit atoms index
        path = [ghost[i] for i in path]
        if not np.all(sum_shift == 0):
            continue
        """
        path = np.array(path)
        i_min = np.argmin(path)
        path = tuple(np.r_[path[i_min:],path[:i_min]])
        # to aline the orientation
        if path[-1] < path[1]:
            path = path[::-1]
            path = tuple(np.r_[path[-1],path[:-1]])
        
        set_rings.add(path)
    
    G.add_edges([(n0, n1)])
    
for i, ring in enumerate(set_rings):
    print(i,ring)

"""
#ring = rings(atoms)

#print(len(ring.rings))
#print(ring.rings[0])

#index = [0, 4, 25, 31, 24, 28, 1, 7]
#for i in index:
#    print(atoms.elements[i], atoms.positions[i])


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