from sovapy.core.file import File
from sovapy.computation.structure_factor import polyhedra
import matplotlib.pyplot as plt
import numpy as np


# Load structural information from a xyz file
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

print("Atom symbols:", atoms.symbols)

# Calculate symmetry measure of polyhedra (q-value),
# whose center is Si and corners are O atoms.
# The threshold distance between the center and a corner 
# to build polyhedora is specified by rmax option.
polys = polyhedra(atoms,center='Si',around='O',rmax=2.0)

# Collect caluculated q-values
q_values = []
for poly in polys:
    if poly.q is not None:
        q_values.append(poly.q)

# Plot the q-value distribution
y, x = np.histogram(q_values, bins=50, range=[0.75,1.01])
plt.bar((x[:-1]+x[1:])*0.5, y, width=0.8*(x[1]-x[0]))
plt.xlim(0.75,1.001)
plt.ylim(0, None)
plt.xlabel('q-values')
plt.ylabel('Counts')
plt.show()

