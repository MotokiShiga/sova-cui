# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 12:04:44 2022

@author: H. Morita
"""

from core.gridding_ctype import Grid
from core.gridding_ctype import metric
import numpy as np
import collections

def tetrahedral(elements,centres,vectors,center,around,r,nth=10):
    MAX_NEIGH = 20
    
    # Set up grid and assign atoms to grid points
    grid = Grid(centres,vectors)
    grid.make_grid()
    
    nmol = centres.shape[0]
    e = collections.Counter(elements).items()
    kind = [k for k, v in e if v > 0]
    ni = [v for k, v in e if v > 0]
    
    _metric = metric(vectors)
    ntypes = len(ni)    
    _rmax = r
    
    #npair = int(ntypes*(ntypes+1)/2)
    
    ncum = []
    ncum.append(ni[0])
    for i in range(1, ntypes):
        ncum.append(ncum[i-1]+ni[i])
    neigh = np.zeros((ntypes,ntypes,MAX_NEIGH), dtype=int)
    
    types = np.zeros(nmol, dtype=int)
    itype = 0
    for i in range(nmol):
        if i == ncum[itype]: itype += 1
        types[i] = itype
        
    centers = np.where(types==center)[0]
    
    c = []
    q = []
    for i in centers:
        grid.neighbours(i, _rmax, 0, nmol-1)
        inei = grid.inei
        nei = len(inei)
        coords = grid.coords
        d = grid.d
        
        n = 0            
        jtype = np.array([types[j] for j in inei])            
        jdx = np.where(jtype == around)[0]
        dd = np.array([d[j] for j in jdx])
        di = np.where(dd <= r)[0]
        
        if len(di) == 4:
            _q = 0.0                
            for j in range(len(di)):
                for k in range(j+1, len(di)):
                    costh = 0.0
                    for ia in range(3):
                        for ib in range(3):
                            costh += _metric[ia][ib]*coords[ia][j]*coords[ib][k]
                    costh = costh/dd[j]/dd[k]
                    costh = max(-1.0, min(costh, 1.0))
                    _q += (costh+1./3.)**2.
            _q = 1.-_q*3./8.
            c.append(i)
            q.append(_q)
            
    return c,q

def calc_tetrahedra(atoms,center,around,r,nth=10):
    MAX_NEIGH = 20
    
    centres = atoms.elementsorted_positions
    vectors = atoms.volume.vectors
    
    # Set up grid and assign atoms to grid points
    grid = Grid(centres,vectors)
    grid.make_grid()
    
    nmol = centres.shape[0]
    #e = collections.Counter(elements).items()
    #kind = [k for k, v in e if v > 0]
    #ni = [v for k, v in e if v > 0]
    ni = atoms.sorted_indexes
    _metric = metric(vectors)
    ntypes = len(ni)    
    _rmax = r
    
    #npair = int(ntypes*(ntypes+1)/2)
    
    print('hhh',center,around)
    
    print(atoms.sorted_elements)
    
    ncum = []
    ncum.append(ni[0])
    for i in range(1, ntypes):
        ncum.append(ncum[i-1]+ni[i])
    neigh = np.zeros((ntypes,ntypes,MAX_NEIGH), dtype=int)
    
    types = np.zeros(nmol, dtype=int)
    itype = 0
    for i in range(nmol):
        if i == ncum[itype]: itype += 1
        types[i] = itype
        
    centers = np.where(types==center)[0]
    
    c = []
    q = []
    for i in centers:
        grid.neighbours(i, _rmax, 0, nmol-1)
        inei = grid.inei
        nei = len(inei)
        coords = grid.coords
        d = grid.d
        
        n = 0            
        jtype = np.array([types[j] for j in inei])            
        jdx = np.where(jtype == around)[0]
        dd = np.array([d[j] for j in jdx])
        di = np.where(dd <= r)[0]
        
        if len(di) == 4:
            _q = 0.0                
            for j in range(len(di)):
                for k in range(j+1, len(di)):
                    costh = 0.0
                    for ia in range(3):
                        for ib in range(3):
                            costh += _metric[ia][ib]*coords[ia][j]*coords[ib][k]
                    costh = costh/dd[j]/dd[k]
                    costh = max(-1.0, min(costh, 1.0))
                    _q += (costh+1./3.)**2.
            _q = 1.-_q*3./8.
            c.append(i)
            q.append(_q)
            
    return c,q