# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 16:58:10 2022

@author: H. Morita
"""

import math
import numpy as np
#for the calculations of the distance matrix and angles (cosine)
from scipy.spatial.distance import pdist, squareform, cosine  

def angle(xi, xj, xk, degree=True):
    """
    calculate angle from 3 vectors / atomic coordinates: i(x,y,z); j(x,y,z); k(x,y,z) 
    xyzarr is the array of all atomic coordinates

    Parameters
    ----------        
    i : list
        neighbor1 position.
    j : list
        center position.
    k : list
        neighbor1 position.

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

def dihedral():
    """
    calculate the dihedral angle from 4 vectors
    atomic coordinates: i(x,y,z); j(x,y,z); k(x,y,z); l(x,y,z)

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    #no warning if division by zero
    np.seterr(invalid='ignore')
    rji = -1*(xyzarr[j] - xyzarr[i])
    rkj = xyzarr[k] - xyzarr[j]
    rlk = xyzarr[l] - xyzarr[k]
    rkj /= np.linalg.norm(rkj)
    v = rji - np.dot(rji, rkj)*rkj
    w = rlk - np.dot(rlk, rkj)*rkj
    x = np.dot(v, w)
    y = np.dot(np.cross(rkj, v), w)
    return np.degrees(np.arctan2(y,x))

def matrix2lattice(m):    
    """
    calculate ftractional to cartesian matrix to 
    lattice parameters (a,b,c,alpha,beta,gamma)

    Parameters
    ----------
    m : numpy matrix
        matrix from ftractional to cartesian.
        
        | a  bcos(gamma) ccos(beta)                                    |
    m = | 0  bsin(gamma) c(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma) |
        | 0  0           cV/sin(gamma)                                 |    
        
        V = sqrt(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2
                 + 2cos(alpha)cos(beta)cos(gamma))
        
    Returns
    -------
    a : float
        a axis length.
    b : float
        b axis length.
    c : float
        c axis length.
    alpha : float
        alpha angle.
    beta : float
        beta angle.
    gamma : TYPE
        gamma angle.

    """
    m = m.T
    a = m[0,0]
    v = m[1,:2]
    b = np.linalg.norm(v)
    cosg = m[1,0]/b
    sing = m[1,1]/b
    v = m[2,:]
    c = np.linalg.norm(v)
    cosb = m[2,0]/c
    cosa = (m[2,0]*cosg+m[2,1]*sing)/c
    alpha = math.acos(cosa)*180./math.pi
    beta = math.acos(cosb)*180./math.pi
    gamma = math.acos(cosg)*180./math.pi    
    
    return a,b,c,alpha,beta,gamma
    