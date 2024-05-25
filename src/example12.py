# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

from core.file import File
from core.molecule import RING

#path = "./sio2_alpha_cristobalite222.cif"
#path = "./data/crystal/sio2_alpha_cristobalite222.cif"
#path = "./data/crystal/sio2_beta_cristobalite.cif"
path = "./data/crystal/sio2_beta_cristobalite222.cif"
path = "./data/crystal/sio2_beta_cristobalite333.cif"
#path = "./data/crystal/sio2_alpha_cristobalite.cif"
f = File.open(path)
atoms = f.getatoms()

#print(atoms.number)

ring = RING()
ring.set_atoms(atoms)

# GUTTMAN, KING, PRIMITIVE
rings = ring.calculate(ring_type=RING.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']],
                       pair_dist_max=[2.0], periodicity=False)

print("Number of non periodic rings : %i" % len(rings))
#for r in ring.rings:
#    print(r)

rings = ring.calculate(ring_type=RING.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']],
                       pair_dist_max=[2.0], periodicity=True)

print('Number of periodic rings : %i' % len(rings))