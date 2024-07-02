# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:26:37 2024

@author: H. Morita
"""

from sova.core.file import File
from sova.computation.rings import RINGs
from sova.computation.cavity import Cavity
from sova.core.data import ResultsFile

path = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()

# Calculate RINGs
ring = RINGs(atoms)
# GUTTMAN, KING, PRIMITIVE
rings = ring.calculate(ring_type=RINGs.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']])

# Calculate Cavity
cavity = Cavity(atoms)
cavity.calculate()

# save file name
path = "./a_SiO2_speed1e11K_rand.hdf5"

# save
with ResultsFile(path, 'w', atoms) as f:
    f.rings = rings
    f.domains = cavity.domains
    f.surface_cavities = cavity.surface_cavities
    f.center_cavities = cavity.center_cavities
    f.flush()

# open
with ResultsFile(path, 'r') as f:
    a = f.atoms
