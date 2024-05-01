# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 12:00:09 2023

@author: H. Morita
"""

import numpy as np
import rmc_cfg
import gridding
import sys

def midpt(path):
    
    # Read starting configuration

    ni, centres, vectors = rmc_cfg.read_rmc_config(path)
    
    n = np.sum(ni)
    ni = ni.tolist()
    ntypes = len(ni)
    ncum = []
    ncum.append(ni[0])
    for i in range(1,ntypes):
        ncum.append(ncum[i-1] + ni[i])
        
    #metric = vectors*vectors
    npar = int(ntypes*(ntypes+1)/2)
    rbond = np.zeros(npar)
    lengths = input(' Maximum bond lengths         : ').split()
    lengths = [float(s) for s in lengths]
    if len(lengths) != npar:
        print(' not match number of lengths {} pair'.format(npar))
    for i, length in enumerate(list(lengths)):
        rbond[i] = length
        
    bonds = input(' Bond               [From,To] : ').split()
    iatom = int(bonds[0])-1
    jatom = int(bonds[1])-1
    patom = int(input(' Type of atom to add          : '))-1
    
    if patom >= ntypes:
        ntypes += 1
        if patom >= ntypes:
            print(' Cannot add that atom type')
            return  
        ni.append(0)
        ncum.append(ncum[-1])
        
    # Find nearest neighbours
    
    grid = gridding.Grid(centres, vectors)
    grid.make_grid()
    
    i1 = min(iatom, jatom)
    j1 = max(iatom, jatom)
    ic = int(i1*(2*ntypes-i1-1)/2+j1)
    itype = 0
    nn = [[] for i in range(n)]
    for i in range(n):
        k = 1
        if i >= ncum[itype]: 
            itype += 1
        if itype != iatom:
            continue
        d, inei = grid.neighbours(i, rbond[ic], 0, n-1)
        for inj in inei:
            jtype = 0
            while inj >= ncum[jtype]:
                jtype += 1
            if jtype == jatom:
                nn[i].append(inj)     

    # Find new atom positions
    
    nadded = 0
    itype = 0
    atoms = centres
    nadd = []
    for i in range(n):
        if i >= ncum[itype]: itype += 1        
        x = atoms[:,i]        
        for j in nn[i]:
            if j <= i: continue
            nadded += 1
            y = atoms[:,j]
            dx = (x + y)/2.
            #if dx[0] == -0.9955831:
            #    print(dx)
            for k in range(3):
                if x[k] < -0.5 and y[k] > 0.5 or y[k] < -0.5 and x[k] > 0.5:
                    dx[k] = dx[k] + 1.
                if dx[k] > 1.: dx[k] = dx[k] - 2.            
            nadd.append(dx)            
    
    # Add new atoms into configuration
    
    atoms = centres.tolist()    
    for i in range(nadded):
        ni[patom] = ni[patom] + 1
        ncum[patom] = ncum[patom] + 1
        for j in range(3): 
            atoms[j].append(0.0)
        n += 1
        
        for j in range(3):
            jtype = ntypes-1            
            shuffling = True
            if patom == jtype: shuffling = False
            while shuffling:                
                atoms[j,ncum[jtype]] = atoms[j,ncum[jtype-1]]
                jtype = jtype-1
                if patom == jtype: 
                    shuffling = False
                        
            atoms[j][ncum[jtype]-1] = nadd[i][j]

    return ni, atoms, vectors

if __name__ == '__main__':
    print(' ADD ATOMS AT BOND CENTRES\n')
    
    ni, centres, vectors = midpt('./si_net.cfg')
    
    path = './si_mid.cfg'
    title = 'title'
    truncated = False
    nmol_types = len(ni)
    nsites = np.ones(nmol_types, dtype=int)
    sites = np.zeros((nmol_types,3))
    neuler = 0
    rmc_cfg.write_rmc_config(path,title,0,0,0,0,
                             truncated,vectors,ni,nsites,sites,
                             neuler,centres)