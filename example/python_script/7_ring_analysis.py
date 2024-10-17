from sovapy.core.file import File
from sovapy.computation.rings import RINGs
import numpy as np
import matplotlib.pyplot as plt

# Load structural information from a cif file
structure_file = "../data/crystal/sio2_beta_cristobalite333.cif"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

# Initialize the class RINGs with structure information
ring = RINGs(atoms)

# Enumerate Guttman rings in the structure
# (Ring types: GUTTMAN, KING, PRIMITIVE)
# Option pair_atom_symbols is not necessary from version 0.5.4.1 for using all chemical bonds.
rings = ring.calculate(ring_type=RINGs.RingType.GUTTMAN)

# rings = ring.calculate(ring_type=RINGs.RingType.GUTTMAN, 
#                        pair_atom_symbols=[['Si', 'O']])

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
"""
Rings with False are closed over the PBD but NOT closed in the real space.
If there are many rings whose outputs are False, the cell size is too small. 
For cif data, the supercell structure can be generated using ase package
using the following code:

from ase.io import read, write

struct = read(structure_file)
struct_new = struct*(2,2,2)  # generate 2x2x2 supercell
write("supercell_structure.cif",struct_new)
""" 

### Statistical analysis  
# extract size, roundness, and roughness of each ring

r_size, r_roundness, r_roughness = list(), list(), list()
for r in rings:
    r_size.append(r.number) # the number of atoms in a ring
    r_roundness.append(r.roundness)
    r_roughness.append(r.roughness)
    
r_size = np.array(r_size)
r_roundness = np.array(r_roundness)
r_roughness = np.array(r_roughness)

#### Ring size distribution
# Maximum ring size
s_max = r_size.max()

# Calculate the histogram of ring size
hist_size = np.zeros(s_max +1, dtype='int')
for s in range(s_max+1):
    hist_size[s] = np.sum(r_size==s)

s_num = np.arange(s_max+1)
plt.figure(figsize=(6,5))
plt.bar(s_num, hist_size)
plt.xlabel('The number of atoms')
plt.ylabel('Counts')
plt.xticks(s_num)
plt.show()


### Roundness and roughness distributions
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.hist(r_roundness, bins=np.linspace(0,1,20))
plt.xlabel('Roundness')
plt.ylabel('Counts')

plt.subplot(1,2,2)
plt.hist(r_roughness, bins=np.linspace(0,1,20))
plt.xlabel('Roughness')
plt.ylabel('Counts')
plt.show()