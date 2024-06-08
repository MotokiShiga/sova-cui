# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from core.file import File
from core.molecule import RING

path = "./data/crystal/sio2_beta_cristobalite.cif"
path = "./data/crystal/sio2_beta_cristobalite333.cif"
f = File.open(path)
atoms = f.getatoms()

ring = RING()
ring.set_atoms(atoms)

# GUTTMAN, KING, PRIMITIVE
rings = ring.calculate(ring_type=RING.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']], 
                       pair_dist_max=[2.0], chain=True)

"""
print('chain type')
for r in rings:
    print(r)

print()
"""  
print('Number of periodic rings : %i' % len(rings))

rings = ring.calculate(ring_type=RING.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']], 
                       pair_dist_max=[2.0])

print('Number of periodic rings : %i' % len(rings))

"""
print('ring type')
for r in rings:
    print(r)
"""