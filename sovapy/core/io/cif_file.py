import numpy as np
from ase import io

from .. import data
from .input_file import InputFile
from .io_utils import FileError

from ...libs.cif2cell.uctools import CellData
import CifFile

class CIFFile(InputFile):
    """ Class for cif file

    Sub-class of `InputFile` to load Crystallographic Information File (CIF) files.

    Attributes
    ----------
    path : string
        Absolute path to a structure file
    cif_grammar : string
        The version of CIF data format
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
        self.cif_grammar = '1.1'
    
    def read_info(self):
        """Read cell data

        Read unit cell (lattice) data from the structure data file

        Raises
        ------
        FileError
            Occur the error when the file cannot be read.
        """
        try:
            self._info.num_frames = 0
            cf = CifFile.ReadCif(self.path, grammar=self.cif_grammar)
            
            # Get the first block
            cb = cf.get(cf.keys()[0])

            # Get cell data
            cell_data = CellData()
            cell_data.HMSymbol = 'P1'
            cell_data.getFromCIF(cb)
            
            # Generate volume data (string) if the crystal system data is available
            self._info.num_frames += 1
            if self._info.num_frames == 1:                
                if cell_data.crystal_system() is not None:
                    # Convert angle data from degree to radian
                    volume_info = 'TRI %f %f %f %f %f %f' % (cell_data.a, cell_data.b, cell_data.c,
                                                              np.radians(cell_data.alpha), 
                                                              np.radians(cell_data.beta),
                                                              np.radians(cell_data.gamma))
                else:
                    # TODO raise exception
                    print('Not found Cell info : cif_file.read_info', cell_data.crystal_system())
                    import sys
                    sys.exit()
                # Save the string of volume information           
                self._info.volumestr = volume_info
            self.has_info = True #Set the flag
        except IOError:
            raise
        except Exception as e:
            raise FileError("Cannot read file info.", e)
    
    def read_atoms(self, frame):
        """Read structure data

        Read structure data from the cif file and generate a Atoms object

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
        FileError
            Error when the file cannot be read.
        """
        try:
            if self.info.num_frames <= frame:
                raise IndexError("Frame {} not found".format(frame))
            
            # Check and read the cif file
            cf = CifFile.ReadCif(self.path, grammar=self.cif_grammar)
            
            # Get blocks
            cfkeys = cf.keys()
            cb = cf.get(cfkeys[0])

            # Get cell data
            cell_data = CellData()
            cell_data.HMSymbol = 'P1'
            cell_data.getFromCIF(cb)
            cell_data.conventional()
            
            # Set cell data got from the file
            if cell_data.HallSymbol != "":
                self.info.volume.crystal_system = cell_data.crystal_system()
                self.info.volume.space_group_number = cell_data.spacegroupnr
                self.info.volume.Hall_symbol = cell_data.HallSymbol
                self.info.volume.Hermann_Mauguin_symbol = cell_data.HMSymbol

            # Load atomic configuration data from the file
            atoms_input = io.read(filename=self.path)            
            symbols = atoms_input.get_chemical_symbols()
            positions = atoms_input.get_positions()
            
            # Set periodic boundary
            self.info.volume.periodic_boundary = [0.0,1.0]

            # Generate Atoms object from the loaded data
            atoms = data.Atoms(positions, None, symbols, self.info.volume)
            # Save structure data described in the structure data
            atoms.original_file_data = data.OriginalStructureData(self.path)

            return atoms
            
        except (IOError, IndexError):
            raise
        except Exception as e:
            raise FileError("Cannot read atom data.", e)       