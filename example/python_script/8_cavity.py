import os
import matplotlib.pyplot as plt
from sovapy.core.file import File
from sovapy.computation.cavity import Cavity

import sovapy
print('sovapy ver: ', sovapy.__version__)

### Load structural information from a xyz file
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K.xyz"
f = File.open(structure_file)

### Get atomic and cell (simulation box) data
atoms = f.get_atoms()


### Cavity analysis
# Initialize Cavity class using atoms
cavity = Cavity(atoms)

# Calculate cavities
# # use default radii
# cavity.calculate(resolution=64)

# Set same radii for all elements
# if gyration_tensor_parameters = True, cavity shape parameters are calculated.
cavity.calculate(resolution=64, cutoff_radii=2.8, gyration_tensor_parameters=True)

# Set custom radii
# cavity.calculate(resolution=64, cutoff_radii={'Si': 2.0, 'O': 3.0}, gyration_tensor_parameters=True)


### Caluculation results
#Three types of cavities: domain, center-based cavity, and surface-based cavity, are calculated. 
#Each of calculated results is obtained by "cavity.domain", "cavity.center_cavities", or "cavity.surface_cavities".
print("The number of domain:")
print(cavity.domains.number)

print("\nThe number of center-based cavities:")
print(cavity.center_cavities.number)

print("\nThe number of surface-based cavities:")
print(cavity.surface_cavities.number)

# Display some properties
print("\n\nCalculate domains")
print('Found {:d} domains'.format(len(cavity.domains.volumes)))
index =  0
print('Domain volume of index {} : {}'.format(index, cavity.domains.volumes[index]))
print('Found {:d} critical domains in file {}. Domain indices: {}'.format(
    len(cavity.domains.critical_domains), os.path.basename(structure_file), 
    cavity.domains.critical_domains))


### Histogram of cavity properties
#
# For the definitions of cavity properties, refer to the "pyMolDyn" paper:  
#  I. Heimbach, F. Rhiem, F. Beule, D. Knodt, J. Heinen, and R.O. Jones.,  
#  "pyMolDyn: Identification, structure, and properties of cavities/vacancies in condensed matter and molecules,"  
#  J. Comput. Chem., 38, 389–394 (2017).   
#   https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.24697

attrs = ['volumes', 'surface_areas', 'squared_gyration_radii', 
         'asphericities', 'acylindricities', 'anisotropies']

plt.figure(figsize=(12,5))
for cnt, attr in enumerate(attrs):
    plt.subplot(2,3,cnt+1)
    vals = getattr(cavity.domains, attr)
    # vals = getattr(cavity.center_cavities, attr)
    # vals = getattr(cavity.surface_cavities, attr)
    plt.hist(vals)
    plt.xlabel(attr)
    plt.ylabel('Counts')
plt.tight_layout()
plt.show()