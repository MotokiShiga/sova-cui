# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 16:58:10 2022

@author: H. Morita
"""

import math
import numpy as np

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
    