# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 16:02:55 2022

@author: H. Morita
"""

import numpy as np
from .. import data
from .input_file import InputFile
from .io_utils import FileError
from ase.io import read
import spglib
from ...libs.cif2cell.uctools import crystal_system
from ...libs.cif2cell.elementdata import ElementData

def get_symbol_from_number(elems, n):
    keys = [k for k, v in elems.items() if v == n]
    if keys:
        return keys[0]
    return None

class EXYZFile(InputFile):
    """
    Implementation on :class:`InputFile` for Open Babel 'exyz' files.
    """
    def __init__(self, path):
        super().__init__(path)
                
        f = open(self.path, "r")
        f.close()
        
    def readinfo(self):
        self._info.num_frames = 0        
        try:
            atoms = read(self.path, format='extxyz')
            dataset = spglib.get_symmetry_dataset(atoms, 
                                                  symprec=1e-5,
                                                  angle_tolerance=-1.0, 
                                                  hall_number=0)
            spacegroup_number = dataset['number']
            _crystal_system = crystal_system(spacegroup_number)
            axes = dataset['std_lattice']
            a, b, c = [np.linalg.norm(axis) for axis in axes ]
            n = np.linalg.norm(axes, axis=0)
            alpha = np.arccos(np.dot(axes[1], axes[2])/n[1]/n[2]) # radian
            beta = np.arccos(np.dot(axes[2], axes[0])/n[2]/n[0])
            gamma = np.arccos(np.dot(axes[0], axes[1])/n[0]/n[1])
                        
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
                print(_crystal_system)
                import sys
                sys.exit()            
            
            self._info.num_frames += 1
            if self._info.num_frames == 1:
                self._info.volumestr = volume_info
                        
            self.inforead = True
        except Exception as e:
            raise FileError("Cannot read file info.", e)
    
    def readatoms(self, frame):
        if self.info.num_frames <= frame:
            raise IndexError("Frame {} not found".format(frame))
        try:
            elems = ElementData()
            symbols = []
            positions = []
            atoms = read(self.path, format='extxyz')
            dataset = spglib.get_symmetry_dataset(atoms, 
                                                  symprec=1e-5,
                                                  angle_tolerance=-1.0, 
                                                  hall_number=0)
            matrix = np.array(self.info.volume.Minv)
            atomic_numbers = dataset['std_types']
            _positions = dataset['std_positions']            
            for i, n in enumerate(atomic_numbers):
                symbol = get_symbol_from_number(elems.elementnr, n)
                pos = np.array(_positions[i])
                pos = matrix.dot(pos)
                symbols.append(symbol)
                positions.append(pos)
            #cx, cy, cz = matrix.dot(np.array([0.5, 0.5, 0.5]))
            #self.info.volume.origin = np.array([-cx, -cy, -cz])            
            self.info.volume.periodic_boundary = [0.0,1.0]
            '''
            for atom in atoms:
                position = (atom.x, atom.y, atom.z)
                symbols.append(atom.symbol)
                positions.append(position)
            '''
            return data.Atoms(positions, None, symbols, self.info.volume)
        except Exception as e:
            raise FileError("Cannot read file info.", e)
    
    