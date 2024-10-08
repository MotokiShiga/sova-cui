from sovapy.core.file import File
from sovapy.computation.rings import RINGs
from sovapy.computation.cavity import Cavity
from sovapy.core.data import ResultsFile


### Load structural information from a xyz file
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(structure_file)


### Get atomic and cell (simulation box) data
atoms = f.getatoms()

### Calculate RINGs
ring = RINGs(atoms)
rings = ring.calculate(ring_type=RINGs.RingType.GUTTMAN, 
                       pair_atom_symbols=[['Si', 'O']])

num_rings1 = len(rings)
s1 = rings[0].number
print("The number of rings: ", num_rings1)
print("The size of the 1st ring:", s1)

# ### Calculate Cavity
cavity = Cavity(atoms)
cavity.calculate()

num_dc1 = cavity.domains.number
v1 = cavity.domains.volumes[0]
print("The number of domain cavities: ", num_dc1)
print("The volume of the 1st domain cavity:", v1)

### Save the structure and computed results 
# Save file name
path = "./a_SiO2_speed1e11K_rand.hdf5"

with ResultsFile(path, 'w', atoms) as f:
    f.rings = rings
    f.domains = cavity.domains
    f.flush()


### Read the file to load calculated results
with ResultsFile(path, 'r') as fr:
    result_atoms = fr.atoms
    result_rings = fr.rings
    result_cavity = fr.domains

num_atoms = len(result_atoms.elements)
print(type(result_atoms))
print("The number of atoms:", num_atoms)

num_rings2 = len(result_rings)
s2 = result_rings[0].number
print(type(result_rings))
print("The number of rings: ", num_rings2)
print("The size of the 1st ring:", s2)

num_dc2 = result_cavity.number
v2 = result_cavity.volumes[0]
print(type(result_cavity))
print("The number of domain cavities: ", num_dc2)
print("The volume of the 1st domain cavity:", v2)




