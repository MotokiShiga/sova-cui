# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 10:01:02 2022

@author: H. Morita
"""

import numpy as np
import ase

from .. import data
from .input_file import InputFile
from .io_utils import FileError

from ...libs.cif2cell.uctools import CellData
from ...libs.cif2cell.elementdata import ElementData
from ...libs.cif2cell.utils import mvmult3

try:
    import CifFile
except ImportError:
    class CifFile:
        """
        Dummy class representing a missing cif2cell module.
        """
        informats = {}

class CIFFile(InputFile):
    """
    Implementation on :class:`InputFile` for Crystallographic Information File (CIF) 'cif' files.
    """
    def __init__(self, path):
        super().__init__(path)
        f = open(self.path, "r")
        f.close()
        
        self.cif_grammar = '1.1'
    
    def readinfo(self):
        try:
            self._info.num_frames = 0
            cf = CifFile.ReadCif(self.path, grammar=self.cif_grammar)
            
            # Get blocks
            cb = cf.get(cf.keys()[0])

            # Get cell data
            cell_data = CellData()
            cell_data.HMSymbol = 'P1'
            cell_data.getFromCIF(cb)
            
            self._info.num_frames += 1
            if self._info.num_frames == 1:                
                # set angle in radian
                if cell_data.crystal_system() is not None:
                    volume_info = 'TRI %f %f %f %f %f %f' % (cell_data.a, cell_data.b, cell_data.c,
                                                              np.radians(cell_data.alpha), 
                                                              np.radians(cell_data.beta),
                                                              np.radians(cell_data.gamma))
                else:
                    # TODO raise exception
                    print('Not found Cell info : cif_file.readinfo', cell_data.crystal_system())
                    import sys
                    sys.exit()
                    
                """
                if cell_data.crystal_system() == 'cubic':                                    
                    volume_info = 'CUB %f' % cell_data.a
                elif cell_data.crystal_system() == 'orthorhombic':
                    volume_info = 'ORT %f %f %f' % (cell_data.a, cell_data.b, cell_data.c)
                elif cell_data.crystal_system() == 'monoclinic':
                    volume_info = 'MON %f %f %f %f' % (cell_data.a, cell_data.b, cell_data.c, 
                                                       np.radians(cell_data.beta))
                elif cell_data.crystal_system() == 'tetragonal':
                    volume_info = 'TET %f %f' % (cell_data.a, cell_data.c)
                elif cell_data.crystal_system() == 'triclinic':
                    volume_info = 'TRI %f %f %f %f %f %f' % (cell_data.a, cell_data.b, cell_data.c,
                                                              np.radians(cell_data.alpha), 
                                                              np.radians(cell_data.beta),
                                                              np.radians(cell_data.gamma))
                elif cell_data.crystal_system() == 'trigonal':
                    volume_info = "RHO %f %f" % (cell_data.a, np.radians(cell_data.alpha))
                elif cell_data.crystal_system() == 'hexagonal':
                    volume_info = 'HEX %f %f' % (cell_data.a, cell_data.c)
                """                
                self._info.volumestr = volume_info
            self.inforead = True
        except IOError:
            raise
        except Exception as e:
            raise FileError("Cannot read file info.", e)
    
    def readatoms(self, frame):         
        try:
            if self.info.num_frames <= frame:
                raise IndexError("Frame {} not found".format(frame))
            
            cf = CifFile.ReadCif(self.path, grammar=self.cif_grammar)
            
            # Get blocks
            cfkeys = cf.keys()
            cb = cf.get(cfkeys[0])

            # Get cell data
            cell_data = CellData()
            cell_data.HMSymbol = 'P1'
            cell_data.getFromCIF(cb)
            #cell_data.primitive()
            cell_data.conventional()
            
            if cell_data.HallSymbol != "":
                self.info.volume.crystal_system = cell_data.crystal_system()
                self.info.volume.space_group_number = cell_data.spacegroupnr
                self.info.volume.Hall_symbol = cell_data.HallSymbol
                self.info.volume.Hermann_Mauguin_symbol = cell_data.HMSymbol                        
            """
            symbols = []
            positions = []
            transmtx = np.array(cell_data.latticevectors)*cell_data.a
            i=0
            for atom in cell_data.atomdata:
                for b in atom:
                    spcsstring = ""
                    for k, v in b.species.items():
                        spcsstring += k+"/"
                    spcsstring = spcsstring.rstrip("/")
                    symbols.append(spcsstring)
                    #v = mvmult3(transmtx.T, b.position)
                    v = transmtx.T.dot(b.position)                    
                    positions.append(np.array([v[0],v[1],v[2]]))
                    i+=1
            """
            atoms_input = ase.io.read(filename=self.path)            
            symbols = atoms_input.get_chemical_symbols()
            positions = atoms_input.get_positions()
            
            self.info.volume.periodic_boundary = [0.0,1.0]
            return data.Atoms(positions, None, symbols, self.info.volume)
            
        except (IOError, IndexError):
            raise
        except Exception as e:
            raise FileError("Cannot read atom data.", e)       