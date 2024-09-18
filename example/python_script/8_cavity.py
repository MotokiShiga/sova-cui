import os
import matplotlib.pyplot as plt
from sovapy.core.file import File
from sovapy.computation.cavity import Cavity

### Load structural information from a xyz file
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(structure_file)

### Get atomic and cell (simulation box) data
atoms = f.getatoms()


### Cavity analysis
# Initialize Cavity class using atoms
cavity = Cavity(atoms)

# Calculate cavities
cavity.calculate()

### Display caluculation results
print("Calculate domains")
print('Found {:d} domains'.format(len(cavity.domains.volumes)))
index =  0
print('Domain volume of index {} : {}'.format(index, cavity.domains.volumes[index]))
print('Found {:d} critical domains in file {}. Domain indices: {}'.format(
    len(cavity.domains.critical_domains), os.path.basename(structure_file), 
    cavity.domains.critical_domains))

plt.figure(figsize=(6,3))
plt.hist(cavity.domains.volumes)
plt.xlabel('Volume')
plt.ylabel('Counts')
plt.show()