from sovapy.core.file import File
from sovapy.computation.rings import RINGs
from sovapy.computation.cavity import Cavity
from sovapy.core.data import ResultsFile

import sovapy
print('sovapy ver: ', sovapy.__version__)

### Load structural information from a xyz file
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K.xyz"
f = File.open(structure_file)

### Get atomic and cell (simulation box) data
atoms = f.getatoms()

num_atoms = len(atoms.elements)
print("The number of atoms:", num_atoms)

### Build chemical bonds
bond_lengths = {('Si', 'O') : 2.0, ('Si', 'Si') : -1, ('O', 'O') : -1}
atoms.set_bond_lengths(bond_lengths)

atoms.bond_summary()
print("")

### Calculate RINGs
ring = RINGs(atoms)
rings_guttman = ring.calculate(ring_type=RINGs.RingType.GUTTMAN)
rings_king = ring.calculate(ring_type=RINGs.RingType.KING)
# rings_primitive = ring.calculate(ring_type=RINGs.RingType.PRIMITIVE)

num_rings1 = len(rings_guttman)
s1 = rings_guttman[0].number
print("")
print("The number of Guttman rings: ", num_rings1)
print("The size of the 1st Guttman ring:", s1)

# ### Calculate Cavity
cavity = Cavity(atoms)
cavity.calculate(resolution=64, cutoff_radii={'Si': 2.0, 'O': 3.0})

num_dc1 = cavity.domains.number
v1 = cavity.domains.volumes[0]
print("")
print("The number of domain cavities: ", num_dc1)
print("The volume of the 1st domain cavity:", v1)

### Save the structure and computed results 
# Save file name
path = "./a_SiO2_speed1e11K_rand.hdf5"

with ResultsFile(path, 'w', atoms=atoms, cavity=cavity) as f:
    f.rings_guttman = rings_guttman
    f.rings_king = rings_king
    # f.rings_primitive = rings_primitive
    f.flush()


### Read the file to load calculated results
with ResultsFile(path, 'r') as fr:
    name                   = fr.name
    version                = fr.version
    result_atoms           = fr.atoms
    result_rings_guttman   = fr.rings_guttman
    result_rings_king      = fr.rings_king
    # result_rings_primitive = fr.rings_primitive
    result_cavity          = fr.cavity

print("\n")
print('Data information:')
print('Package: ', name)
print('Version: ', version)

num_atoms = len(result_atoms.elements)
print("")
print(type(result_atoms))
print("The number of atoms:", num_atoms)
print("Atom symbols:", result_atoms.symbols)
print("")
result_atoms.bond_summary()

num_rings2 = len(result_rings_guttman)
s2 = result_rings_guttman[0].number
print("")
print(type(result_rings_guttman))
print("The number of Guttman rings: ", num_rings2)
print("The size of the 1st Guttman ring:", s2)

num_dc2 = result_cavity.domains.number
v2 = result_cavity.domains.volumes[0]
print(type(result_cavity))
print("")
print("The number of domain cavities: ", num_dc2)
print("The volume of the 1st domain cavity:", v2)
print("Calculation settings: resolution={:}, cutoff_radii={:}".format(result_cavity.resolution,result_cavity.cutoff_radii))




