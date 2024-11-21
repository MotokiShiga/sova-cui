from sovapy.core.file import File
from sovapy.core.analysis import (PDFAnalysis, BondAngleAnalysis, CoordinationNumberAnalysis,
                                  TetrahedralOrderAnalysis, RingAnalysis, CavityAnalysis)


import sovapy
print('sovapy ver: ', sovapy.__version__)

# Load structural information from a xyz file
# The second line (CUB 24.713) in the xyz file indicates the shape of cell and its length.
# (CUB means cubics.)
structure_file = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()
print("Atom symbols:", atoms.symbols)

# Setting of maximum bond lengths for atomic element pairs
# Set -1 to pairs for which you don't want to build bonds.
print('Set bond length:')
bond_lengths = {('Si', 'O') : 2.0, ('Si', 'Si') : -1, ('O', 'O') : -1}
atoms.set_bond_lengths(bond_lengths)
print('\n\n')

# PDF analysis
print('PDF (Pair Distribution Function) Analysis:')
pdf = PDFAnalysis(atoms, dr=0.05, dq=0.05, qmin=0.3, qmax=25.0)
pdf.run()
pdf.plot()
pdf.save_atoms_to_hdf5(mode="w") # the first analysis should save Atoms object using the method
pdf.save_to_hdf5()
print('\n\n')

# Coordination Number Analysis
print('Coordination Number Analysis')
coord_num = CoordinationNumberAnalysis(atoms)
coord_num.run()
coord_num.plot()
coord_num.save_to_hdf5()
print('\n\n')

# Bond angle analysis
print('Bond angle analysis:')
bond_angle = BondAngleAnalysis(atoms, bins=50)
bond_angle.run()
bond_angle.plot()
bond_angle.save_to_hdf5()
print('\n\n')

# Tetrahedral order analysis
print('Tetrahedral order analysis:')
list_cc_dist = [['Si','O',2.0],['Si','Si',3.5]]
tetra = TetrahedralOrderAnalysis(atoms, bins=100, list_cc_dist=list_cc_dist)
tetra.run()
tetra.plot()
tetra.save_to_hdf5()
print('\n\n')

# Ring analysis
print('Ring analysis:')
ring = RingAnalysis(atoms, guttman=True, king=True, primitive=False, cutoff_primitive=24)
ring.run()
ring.plot()
ring.save_to_hdf5()
print('\n\n')

# Cavity analysis
print('Cavity analysis:')
# cutoff_radii = {'Si': 2.8, 'O': 2.9} # radius setting for each atoms
cutoff_radii = 2.8 # same radii for all atoms
cavity = CavityAnalysis(atoms, resolution=64, cutoff_radii=cutoff_radii)
cavity.run()
cavity.plot()
cavity.save_to_hdf5()
print('\n\n')
