#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 10:27:36 2023

@author: morita
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform, cosine    #for the calculations of the distance matrix and angles (cosine)

def calc_angle(xi, xj, xk, degree=True):
    """
    calculate angle from 3 vectors / atomic coordinates: i(x,y,z); j(x,y,z); k(x,y,z) 
    xyzarr is the array of all atomic coordinates

    Parameters
    ----------
    xyzarr : list
        positions list.
    i : int
        index.
    j : int
        index.
    k : int
        index.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    rij = xi - xj 
    rkj = xk - xj
	#remove if cosine fails
	#cos_theta = np.dot(rij, rkj)
	#sin_theta = np.linalg.norm(np.cross(rij, rkj))
	#theta = np.arctan2(sin_theta, cos_theta)
	#scipy pdist cosine instead of the 3 lines above 
    theta = cosine(rij, rkj)
    if degree == False:
        return 1.-theta
    theta = np.arccos(1.-theta) 
    return np.degrees(theta)

class Polyhedra(object):
    def __init__(self, center=None, arounds=None):
        self.center = center
        self.arounds = arounds
        self.q = 0.
        
    def calc_order(self):
        self.q = 0.0 
        for j, around in enumerate(self.arounds) :
            for k in range(j+1, len(self.arounds)):
                costh = calc_angle(around, self.center, self.arounds[k], degree=False)
                self.q += (costh+1./3.)**2.                
        self.q = 1. - self.q*3./8.
