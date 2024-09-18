#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 10:27:36 2023

@author: morita
"""

import numpy as np
from .utils import angle
from scipy.spatial import Delaunay

class Polyhedron(object):
    def __init__(self, atoms,center,neighbors,vectors,shifts):
        self.atoms = atoms
        self.center = center
        self.neighbors = neighbors
        self.vectors = vectors # "relative vector from center position
        self.shifts = shifts
        self._volume = None
        self._q = None        
    
    @property
    def q(self):
        """
        calculate tetrahedral order only tetrahedron.

        Returns
        -------
        float 
            tetrahedral order.

        """
        if self._q is None and len(self.neighbors) == 4:            
            self._q = 0.0            
            position = np.zeros(3)
            for j, neighbor in enumerate(self.neighbors):            
                for k in range(j+1, len(self.neighbors)):
                    costh = angle(self.vectors[j],position, self.vectors[k], degree=False)
                    self._q += (costh+1./3.)**2.         
            self._q = 1. - self._q*3./8.
        return self._q
    
    @property
    def volume(self):
        if self._volume is None:
            self._volume = 0.
            if len(self.vectors) > 3:
                points = np.array(self.vectors)
                tetrahedra = Delaunay(points).simplices
                for t in tetrahedra:
                    v0 = self.vectors[1]-self.vectors[0]
                    v1 = self.vectors[2]-self.vectors[0]
                    v2 = self.vectors[3]-self.vectors[0]                
                    cross_product = np.cross(v1-v0,v2-v0)
                    v = np.linalg.norm(cross_product)/6
                    self._volume += v
        return self._volume