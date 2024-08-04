# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita and M. Shiga
"""

from sova.core.file import File
from sova.computation.rings import RINGs

# Load structural information from a cif file
structure_file = "./data/crystal/sio2_beta_cristobalite333.cif"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

# Initialize the class RINGs with structure information
ring = RINGs(atoms)

# Enumerate Guttman rings in the structure
# (Ring types: GUTTMAN, KING, PRIMITIVE)
rings = ring.calculate(ring_type=RINGs.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']])

"""
Parallel computation is available using num_parallel option as follows:
rings = ring.calculate(ring_type=RINGs.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']], num_parallel=8)

To use the maximum numbers of CPU cores, set -1.
"""

# Computed result
print("The number of rings in the structure : %i" % len(rings))

n = 2 # Ring ID to output computated result
print("The number of atoms in the n-th ring:")
print(rings[n].number)

print("Roundness and roughness in the n-th ring:")
print([rings[n].roundness, rings[2].roughness])

print("Does the ring cross the boundary of the cell?")
print(rings[n].over_boundary)

print("Does the ring is closed in the real space?")
print(rings[n].close)
# If there are many rings whose outputs are True, the cell size is too small. 
# For cif data, the supercell structure can be generated using ase package
# using the following code:
#
"""
from ase.io import read, write

struct = read(structure_file)
struct_new = struct*(2,2,2)  # generate 2x2x2 supercell
write("supercell_structure.cif",struct_new)
""" 

