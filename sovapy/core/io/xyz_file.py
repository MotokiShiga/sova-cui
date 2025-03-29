import numpy as np
from .. import data
from .input_file import InputFile
from .io_utils import FileError
from ..volumes import HexagonalVolume,NonVolume

class XYZFile(InputFile):
    """Class for XYZ files

    Sub-class of `InputFile` to load xyz files.
    
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
        try:
            self._info.num_frames = 0
            with open(self.path.encode("utf-8"), 'r') as f:
                try:
                    num_atoms = f.readline().replace("\n", "")
                    if not num_atoms:
                        print("Warning: Not found the number of atoms in the xyz file!")
                        raise Exception
                    num_atoms = int(num_atoms)
                    #Read cell info (volume)
                    volume_info = f.readline().replace("\n", "")
                    if not volume_info:
                        print("Warning: Not found cell information in the xyz file!")
                        raise Exception
                    # Read the following lines to check the file
                    for i in range(num_atoms):
                        f.readline()
                    self._info.num_frames += 1
                    #Save cell info (volume) if the data is valid
                    if self._info.num_frames == 1:
                        self._info.volumestr = volume_info
                except StopIteration:
                    pass
            # Update the flag
            self.has_info = True
        except IOError:
            raise
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
            Error related with frame index specification
        StopIteration
            _description_
        IndexError
            _description_
        """
        try:
            if self.info.num_frames <= frame:
                raise IndexError("Frame {} not found".format(frame))
            with open(self.path.encode("utf-8"), 'r') as f:
                try:
                    # Skip the first frames
                    for i in range(frame):
                        num_atoms = f.readline().replace("\n", "")
                        if not num_atoms:
                            print("Warning: Not found the number of atoms in the xyz file!")
                            raise Exception
                        num_atoms = int(num_atoms)
                        f.readline()
                        for i in range(num_atoms):
                            f.readline()

                    # Read atomic configuration in the specified frame
                    symbols, positions = [], []
                    
                    # Read the number of atoms
                    num_atoms = f.readline().replace("\n", "")
                    if not num_atoms:
                        raise StopIteration
                    num_atoms = int(num_atoms)

                    # Load volumestr but not saved here 
                    # because is saved in the read_info method
                    volume_info = f.readline() # volumestr
                    if not volume_info:
                        print("Warning: Not found cell information in the xyz file!")
                        raise Exception

                    # Read atomic symbols and positions
                    for i in range(num_atoms):
                        line = f.readline()
                        if line.strip():
                            symbol, x, y, z = line.split()[:4]
                            position = (float(x), float(y), float(z))
                            symbols.append(symbol)
                            positions.append(position)

                    # Set the origin of the coordinate
                    if isinstance(self.info.volume, HexagonalVolume):
                        cx,cy,cz = 0.0, 0.0, 0.0
                    elif isinstance(self.info.volume, NonVolume):
                        cx,cy,cz = 0.5*(np.max(positions,axis=0)+np.min(positions,axis=0))
                    else:
                        cx,cy,cz = np.array(self.info.volume.Minv).dot(np.array([0.5,0.5,0.5]))
                        self.info.volume.origin = np.array([-cx,-cy,-cz])

                    # Generate Atoms object from the loaded data
                    atoms = data.Atoms(positions, None, symbols, self.info.volume)
                    # Save structure data described in the structure data
                    atoms.original_file_data = data.OriginalStructureData(self.path)
                except StopIteration:
                    raise IndexError("Frame {} not found".format(frame))
            return atoms
        except (IOError, IndexError):
            raise
        except Exception as e:
           raise FileError("Cannot read atom data.", e)
            
    #Deprecated
    # @staticmethod
    # def write(result,filepath):
    #     with open(filepath, mode='w') as f:
    #         f.write("{:<10d}\n".format(result.atoms.number))
    #         f.write("\n")
    #         for i in range(result.atoms.number):
    #             f.write('{:<2s}'.format(result.atoms.elements[i]))
    #             f.write('{:12.07f}'.format(result.atoms.positions[i][0]))
    #             f.write('{:12.07f}'.format(result.atoms.positions[i][1]))
    #             f.write('{:12.07f}\n'.format(result.atoms.positions[i][2]))
            