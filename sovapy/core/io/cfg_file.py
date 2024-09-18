# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 21:10:49 2022

@author: H. Morita
"""

import os
import numpy as np
from .. import data
from .input_file import InputFile
from .io_utils import FileError
from ...computation.utils import matrix2lattice
import collections

class CFGFile(InputFile):
    """
    Implementation on :class:`InputFile` for RMC 'cfg' files.
    """
    def __init__(self, path):
        super().__init__(path)
        
        self.elements = None
        
        f = open(self.path, "r")
        f.close()
    
    def readinfo(self):
        try:
            self._info.num_frames = 0
            with open(self.path.encode("utf-8"), 'r') as f:
                try:
                    while True:
                        ver = f.readline()
                        title = f.readline()

                        line = f.readline()
                        item = line.split()
                        while len(item) == 0:
                            line = f.readline()
                            item = line.split()
                                                
                        ngen, ntried, nacc = item[:3]
                        ngen = int(ngen)
                        ntried = int(ntried)
                        nacc = int(nacc)
                        nsaved = int(f.readline().split()[0])
                        
                        f.readline()
                        num_atoms = int(f.readline().split()[0])
                        nmol_types = int(f.readline().split()[0])
                        nsites_max = int(f.readline().split()[0])
                        neuler = int(f.readline().split()[0])
                        
                        f.readline()
                        items = f.readline().split()
                        if items[0] == 'F':
                            truncated = False
                        elif items[0] == 'T':
                            truncated = True
                        f.readline()
                        
                        vectors = np.identity(3)
                        for i in range(3):
                            line = f.readline()
                            xyz = np.fromstring(line, dtype=float, sep=' ') 
                            vectors[i] = xyz

                        ni = []
                        nsites = []
                        sites = []
                        f.readline()
                        for i in range(nmol_types):
                            line = f.readline()
                            ni.append(int(line.split()[0]))
                            line = f.readline()
                            n = int(line.split()[0])
                            nsites.append(n)
                            s = []
                            for k in range(n):
                                line = f.readline()
                                s.append(np.fromstring(line, dtype=float, sep=' '))
                            sites.append(s)
                            f.readline()
                        
                        # TODO information add 
                        self._info.nmol_types = nmol_types
                        self._info.ni = ni

                        self._info.num_frames += 1
                        if self._info.num_frames == 1:                            
                            m = vectors*2
                            a, b, c, alpha, beta, gamma = matrix2lattice(m)
                            angle = np.array([alpha, beta, gamma])
                            
                            volume_info = None
                            if a == b and b == c and  alpha == 90. and alpha == beta and beta == gamma: # Cubic
                                volume_info = 'CUB %f' % a
                            elif a != b and b != c and c != a and alpha == 90. and alpha == beta and beta == gamma: # Orthorhombic
                                volume_info = 'ORT %f %f %f' % (a, b, c)
                            elif a != b and b != c and c != a and alpha == 90. and alpha != beta and alpha == gamma: # Monoclinic
                                volume_info = 'MON %f %f %f %f' % (a, b, c, np.radians(beta))
                            elif a != b and a == c and alpha == 90. and alpha == beta and beta == gamma: # Tetragonal
                                volume_info = 'TET %f %f' % (a, c)
                            elif a != b and b != c and c != a and alpha != beta and beta != gamma and gamma != alpha: # Triclinic
                                volume_info = 'TRI %f %f %f %f %f %f' % (a, b, c, np.radians(alpha), np.radians(beta), np.radians(gamma))
                            elif a == b and b == c and alpha != 90. and alpha == beta and beta == gamma: # Trigonal
                                volume_info = "RHO %f %f" % (a, np.radians(alpha))
                            else:
                                pass
                            
                            self._info.volumestr = volume_info
                        break

                except StopIteration:
                    pass
            self.inforead = True
        except IOError:
            raise
        except Exception as e:
            raise FileError("Cannot read file info.", e)
            
    def readelements(self,nmol_types):
        elem_file = os.path.splitext(self.path)[0] + '.elm'
        self.elements = None
        if os.path.exists(elem_file) == True:
            try:
                self.elements = []
                with open(elem_file.encode("utf-8"), 'r') as f:
                    for i in range(nmol_types):
                        self.elements.append(f.readline().strip())
            except FileNotFoundError:
                print('File Not Found : ', elem_file)
    
    def readatoms(self, frame,elements=None):
        self.elements = elements 
        try:
            if self.info.num_frames <= frame:
                raise IndexError("Frame {} not found".format(frame))
            with open(self.path.encode("utf-8"), 'r') as f:
                try:
                    # Skip the first frames
                    for i in range(frame):
                        num_atoms = f.readline()

                    # actually read the molecule
                    symbols = []
                    positions = []
                    ver = f.readline()
                    title = f.readline()

                    line = f.readline()
                    item = line.split()
                    while len(item) == 0:
                        line = f.readline()
                        item = line.split()
                    
                    ngen, ntried, nacc = item[:3]
                    ngen = int(ngen)
                    ntried = int(ntried)
                    nacc = int(nacc)
                    nsaved = int(f.readline().split()[0])
                    
                    f.readline()
                    num_atoms = int(f.readline().split()[0])
                    nmol_types = int(f.readline().split()[0])
                    nsites_max = int(f.readline().split()[0])
                    neuler = int(f.readline().split()[0])
                    
                    f.readline()
                    item = f.readline().split()[0]
                    if item == 'F':
                        truncated = False
                    elif item == 'T':
                        truncated = True
                    f.readline()
                    
                    vectors = np.identity(3)
                    self.info.volume._vectors = vectors
                    for i in range(3):
                        line = f.readline()
                        xyz = np.fromstring(line, dtype=float, sep=' ') 
                        vectors[i] = xyz
                    
                    ni = []
                    nsites = []
                    sites = []
                    f.readline()
                    for i in range(nmol_types):
                        line = f.readline()
                        ni.append(int(line.split()[0]))
                        line = f.readline()
                        n = int(line.split()[0])
                        nsites.append(n)
                        s = []
                        for k in range(n):
                            line = f.readline()
                            s.append(np.fromstring(line, dtype=float, sep=' '))
                        sites.append(s)
                        f.readline()
                    
                    #self.readelements(nmol_types)
                                    
                    # positions
                    for i in range(nmol_types):
                        for j in range(ni[i]):
                            symbol = self.elements[i]
                            symbols.append(symbol)
                            line = f.readline()
                            x, y, z = line.split()
                            v = np.array([float(x), float(y), float(z)])
                            position = np.dot(v, vectors)
                            positions.append(position)
                    
                    # set origin point
                    self.info.volume.origin[0] = -vectors[0][0]
                    self.info.volume.origin[1] = -vectors[1][1]
                    self.info.volume.origin[2] = -vectors[2][2]
                    
                except StopIteration:
                    raise IndexError("Frame {} not found".format(frame))
            return data.Atoms(positions, None, symbols, self.info.volume, True)
        except (IOError, IndexError):
            raise
        except Exception as e:
            raise FileError("Cannot read atom data.", e)
            
    @staticmethod
    def write(result,filepath):
        with open(filepath, mode='w') as f:
            f.write("(Version 3 format configuration file)\n")
            f.write(" From " + result.filepath)
            f.write('\n\n')            
            f.write(" ")
            f.write("{:>10d}".format(0))
            f.write("{:>10d}".format(0))
            f.write("{:>10d}".format(0))
            f.write(" moves generated, tried, accepted\n")
            f.write(" ")
            f.write("{:>10d}".format(0))
            f.write("                     configurations saved\n\n")
            f.write(" ")
            nmol = result.atoms.number
            vectors = result.atoms.volume.vectors
            ni = []            
            for k, v in collections.Counter(result.atoms.elements).items():
                ni.append(v)
            nmol_types = len(ni)
            truncated = False
            #print(c)            
            f.write("{:>10d}".format(nmol))
            f.write(" molecules (of all types)\n")
            f.write(" ")
            f.write("{:>10d}".format(nmol_types))
            f.write(" types of molecule\n")
            f.write(" ")
            f.write("{:>10d}".format(0))
            f.write(" is the largest number of atoms in a molecule\n")
            f.write(" ")
            f.write("{:>10d}".format(0))
            f.write(" Euler angles are provided\n\n")
            f.write("          ")
            if truncated == True:
                f.write('T')
            else:
                f.write('F')
            f.write(" (Box is not truncated octahedral)\n")
            f.write("           ")
            f.write(" Defining vectors are:\n")            
            for i in range(3):
                f.write("           ")
                f.write('{:11.06f}'.format(round(vectors[i][0],7)))
                f.write('{:11.06f}'.format(round(vectors[i][1],7)))
                f.write('{:11.06f}\n'.format(round(vectors[i][2],7)))
            f.write('\n')
            
            for i in range(nmol_types):
                f.write("{:>11d}".format(ni[i]))
                f.write(" molecules of type")
                f.write("{:>3d}\n".format(i+1))
                f.write("{:>11d}".format(1))
                f.write(" atomic sites\n")
                f.write("           ")
                f.write('{:11.06f}'.format(0))
                f.write('{:11.06f}'.format(0))
                f.write('{:11.06f}'.format(0))
                f.write('\n\n')
            
            for i in range(nmol):
                f.write('{:11.07f}'.format(result.atoms.norm_positions[i][0]))
                f.write('{:11.07f}'.format(result.atoms.norm_positions[i][1]))
                f.write('{:11.07f}\n'.format(result.atoms.norm_positions[i][2]))