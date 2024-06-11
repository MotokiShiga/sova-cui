# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname('__file__'), '..'))

from core.file import File
from core.molecule import RINGs

path = "../data/crystal/sio2_beta_cristobalite333.cif"
f = File.open(path)
atoms = f.getatoms()

ring = RINGs(atoms)

# GUTTMAN, KING, PRIMITIVE
rings = ring.calculate(ring_type=RINGs.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']])

print("Number of rings : %i" % len(rings))

print(rings[2].close)
print(rings[2].roundness)
print(rings[2].roughness)
