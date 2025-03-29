import ase, spglib
import numpy as np

from .. import data
from .input_file import InputFile
from .io_utils import FileError
from ...libs.cif2cell.uctools import crystal_system
from ...libs.cif2cell.elementdata import ElementData

def get_symbol_from_number(elems, n):
    """Convert atomic numbers to symbols

    Parameters
    ----------
    elems : _type_
        _description_
    n : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    keys = [k for k, v in elems.items() if v == n]
    if keys:
        return keys[0]
    return None

class EXYZFile(InputFile):
    """Class for Extended XYZ files

    Sub-class of `InputFile` to load extended xyz (exyz or extxyz) files.
    
    Attributes
    ----------
    path : string
        Absolute path to a structure file
    _info : core.data.FileInfo object
        Meta data including #frames, unit cell info (volume)
    has_info : bool
        If the structure infomation has been already read
    """
    def __init__(self, path):
        """Initialization

        Parameters
        ----------
        path : string
            Absolute path to a structure file
        """
        super().__init__(path)
        f = open(self.path, "r")
        f.close()
        
    def read_info(self):
        """Read cell data

        Read unit cell (lattice) data from the structure data file

        Raises
        ------
        FileError
            Error when the file cannot be read.
        """
        self._info.num_frames = 0        
        try:
            # Load the structure data
            atoms = ase.io.read(self.path, format='extxyz')
            # Find symmetry operations from the structure
            dataset = spglib.get_symmetry_dataset(atoms, 
                                                  symprec=1e-5,
                                                  angle_tolerance=-1.0, 
                                                  hall_number=0)
            spacegroup_number = dataset['number']
            _crystal_system = crystal_system(spacegroup_number)
            # Lattice vectors
            axes = dataset['std_lattice']
            # Compute Lattice constants (three lengths)
            a, b, c = [np.linalg.norm(axis) for axis in axes ]
            # Compute Lattice constants (three angles (radian))
            n = np.linalg.norm(axes, axis=0)
            alpha = np.arccos(np.dot(axes[1], axes[2])/n[1]/n[2])
            beta = np.arccos(np.dot(axes[2], axes[0])/n[2]/n[0])
            gamma = np.arccos(np.dot(axes[0], axes[1])/n[0]/n[1])

            # Generate a string about the volume (cell) data       
            if _crystal_system == 'cubic':                                    
                volume_info = 'CUB %f' % a
            elif _crystal_system == 'orthorhombic':
                volume_info = 'ORT %f %f %f' % (a, b, c)
            elif _crystal_system == 'monoclinic':
                volume_info = 'MON %f %f %f %f' % (a, b, c, beta)
            elif _crystal_system == 'tetragonal':
                volume_info = 'TET %f %f' % (a, c)
            elif _crystal_system == 'triclinic':
                volume_info = 'TRI %f %f %f %f %f %f' % (a, b, c, alpha, beta, gamma)
            else:
                print('Error: Crystal_system was identifyed as ', _crystal_system)
                raise Exception
            
            self._info.num_frames += 1
            if self._info.num_frames == 1:
                # Set volume (cell) data
                self._info.volumestr = volume_info
            # Update the flag           
            self.has_info = True
        except Exception as e:
            raise FileError("Cannot read file info.", e)
    
    def read_atoms(self, frame):
        """Read structure data

        Read structure data from the xyz file and generate a Atoms object

        Parameters
        ----------
        frame : int
            Frame index

        Returns
        -------
        atoms : core.data.Atoms object
            Atoms object that contains cell and atomic configuration data

        Raises
        ------
        IndexError
            _description_
        FileError
            _description_
        """
        if self.info.num_frames <= frame:
            raise IndexError("Frame {} not found".format(frame))
        try:
            # Load structure data using ase package
            atoms = ase.io.read(self.path, format='extxyz')
            # Find symmetry operations from the structure
            dataset = spglib.get_symmetry_dataset(atoms, 
                                                  symprec=1e-5,
                                                  angle_tolerance=-1.0, 
                                                  hall_number=0)
            # The fracional-to-cartesian transformation matrix
            matrix_F2C = np.array(self.info.volume.Minv)
            # Get atomic numbers and fractional coordinates (relative positions)
            atomic_numbers = dataset['std_types']
            _positions = dataset['std_positions']
            # Convert the data to atomic symbols and cartesian coodinates
            elems = ElementData()
            symbols, positions = [], []
            for i, n in enumerate(atomic_numbers):
                symbol = get_symbol_from_number(elems.elementnr, n)
                pos = np.array(_positions[i])
                pos = matrix_F2C.dot(pos) #convert from fractional coodinate to cartesian coodinate
                symbols.append(symbol)
                positions.append(pos)
            self.info.volume.periodic_boundary = [0.0,1.0]

            # Generate Atoms object from the stcuture data
            atoms = data.Atoms(positions, None, symbols, self.info.volume)
            # Save structure data described in the structure data
            atoms.original_file_data = data.OriginalStructureData(self.path)
            return atoms
        except Exception as e:
            raise FileError("Cannot read file info.", e)
    
    