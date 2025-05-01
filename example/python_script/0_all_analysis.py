from sovapy.core.file import File
from sovapy.core.analysis import (PDFAnalysis, BondAngleAnalysis, CoordinationNumberAnalysis,
                                  TetrahedralOrderAnalysis, RingAnalysis, CavityAnalysis)
from sovapy.core.data import ResultsFile

import multiprocessing

import sovapy

def main():
    print('sovapy ver: ', sovapy.__version__)

    # Load structural information from a xyz file
    # The second line (CUB 24.713) in the xyz file indicates the shape of cell and its length.
    # (CUB means cubics.)
    structure_file = "../data/amorphous_md/a_SiO2_speed1e11K.xyz"
    f = File.open(structure_file)

    # Get atomic and cell (simulation box) data
    atoms = f.get_atoms()
    print("Atom symbols:", atoms.symbol_set)

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
    tetra = TetrahedralOrderAnalysis(atoms, num_bins=100, list_cc_dist=list_cc_dist)
    tetra.run()
    tetra.plot()
    tetra.save_to_hdf5()
    print('\n\n')

    # Ring analysis
    print('Ring analysis:')
    # To use only single cpu core 
    ring = RingAnalysis(atoms, guttman=True, king=True, primitive=False, cutoff_primitive=24)
    # To use multiple cpu cores 
    # ring = RingAnalysis(atoms, guttman=True, king=True, primitive=False, cutoff_primitive=24, num_parallel=-1)
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

    # Load Computation results
    path = './a_SiO2_speed1e11K.hdf5'
    with ResultsFile(path, 'r') as fr:
        name                   = fr.name
        version                = fr.version
        result_atoms           = fr.atoms
        result_rings_guttman   = fr.rings_guttman
        result_rings_king      = fr.rings_king
        # result_rings_primitive = fr.rings_primitive
        result_cavity          = fr.cavity

    print('Data information:')
    print('Package: ', name)
    print('Version: ', version)

    num_rings2 = len(result_rings_guttman)
    s2 = result_rings_guttman[0].size
    print(type(result_rings_guttman))
    print("The number of Guttman rings: ", num_rings2)
    print("The size of the 1st Guttman ring:", s2)

    num_dc2 = result_cavity.domains.number
    v2 = result_cavity.domains.volumes[0]
    print(type(result_cavity))
    print("The number of domain cavities: ", num_dc2)
    print("The volume of the 1st domain cavity:", v2)
    print("Calculation settings: resolution={:}, cutoff_radii={:}".format(result_cavity.resolution,result_cavity.cutoff_radii))


if __name__ == '__main__':
    # Necessary for parallel computation
    multiprocessing.freeze_support()
    
    main()