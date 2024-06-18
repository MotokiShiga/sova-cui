# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 16:21:46 2022

@author: H. Morita
"""
import numpy as np
from .. import data
from .input_file import InputFile
from .io_utils import FileError
from ..volumes import HexagonalVolume,NonVolume

class XYZFile(InputFile):
    """
    Implementation on :class:`InputFile` for Open Babel 'xyz' files.
    """
    def __init__(self, path):
        super().__init__(path)
        f = open(self.path, "r")
        f.close()
    
    def readinfo(self):
        try:
            self._info.num_frames = 0
            with open(self.path.encode("utf-8"), 'r') as f:
                try:
                    while True:
                        num_atoms = f.readline().replace("\n", "")
                        if not num_atoms:
                            break
                        num_atoms = int(num_atoms)
                        volume_info = f.readline().replace("\n", "")
                        for i in range(num_atoms):
                            f.readline()
                        self._info.num_frames += 1
                        if self._info.num_frames == 1:
                            self._info.volumestr = volume_info
                except StopIteration:
                    pass
            self.inforead = True
        except IOError:
            raise
        except Exception as e:
            raise FileError("Cannot read file info.", e)        

    def readatoms(self, frame):
        try:
            if self.info.num_frames <= frame:
                raise IndexError("Frame {} not found".format(frame))
            with open(self.path.encode("utf-8"), 'r') as f:
                try:
                    # Skip the first frames
                    for i in range(frame):
                        num_atoms = f.readline().replace("\n", "")
                        if not num_atoms:
                            break
                        num_atoms = int(num_atoms)
                        f.readline()
                        for i in range(num_atoms):
                            f.readline()
                    # actually read the molecule
                    symbols = []
                    positions = []
                    num_atoms = f.readline().replace("\n", "")
                    if not num_atoms:
                        raise StopIteration
                    num_atoms = int(num_atoms)
                    f.readline() # volumestr
                    for i in range(num_atoms):
                        line = f.readline()
                        if line.strip():
                            symbol, x, y, z = line.split()[:4]
                            position = (float(x), float(y), float(z))
                            symbols.append(symbol)
                            positions.append(position)
                    # origin
                    if isinstance(self.info.volume, HexagonalVolume):
                        cx,cy,cz = 0.0, 0.0, 0.0
                    elif isinstance(self.info.volume, NonVolume):
                        cx,cy,cz = 0.5*(np.max(positions,axis=0)+np.min(positions,axis=0))
                    else:
                        cx,cy,cz = np.array(self.info.volume.Minv).dot(np.array([0.5,0.5,0.5]))
                        self.info.volume.origin = np.array([-cx,-cy,-cz])
                except StopIteration:
                    raise IndexError("Frame {} not found".format(frame))
            return data.Atoms(positions, None, symbols, self.info.volume)
        except (IOError, IndexError):
            raise
        #except Exception as e:
        #    raise FileError("Cannot read atom data.", e)
            
    @staticmethod
    def write(result,filepath):
        with open(filepath, mode='w') as f:
            f.write("{:<10d}\n".format(result.atoms.number))
            f.write("\n")
            for i in range(result.atoms.number):
                f.write('{:<2s}'.format(result.atoms.elements[i]))
                f.write('{:12.07f}'.format(result.atoms.positions[i][0]))
                f.write('{:12.07f}'.format(result.atoms.positions[i][1]))
                f.write('{:12.07f}\n'.format(result.atoms.positions[i][2]))
            