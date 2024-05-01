# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 09:23:10 2023

@author: H. Morita
"""

import os
import sys
import numpy as np
import rmc_cfg 
from gridding import Grid
from distutils.util import strtobool 

def sign(a, b):
    if(b >= 0.0):
        return np.fabs(a)
    else:
        return -np.fabs(a)

def moveout(ni, centres, vectors):
    
    max_atoms = 30000
        
    truncated = False
    ntypes = len(ni)
    n = np.sum(ni)
    ind = np.zeros((max_atoms, ntypes), dtype=int)    
    ncum = np.zeros(ntypes)
    ncum[0] = ni[0]
    for i in range(1,ntypes):
        ncum[i] = ncum[i-1] + ni[i]
    
    np.random.seed(12345)
    metric = vectors*vectors
    '''
    metric = np.identity(3)
    for j in range(3):
        for i in range(3):
            metric[i,j] = 0.
            for k in range(3):
                metric[i,j] += vectors[k,i]*vectors[k,j]
    '''
    npar = int(ntypes*(ntypes+1)/2)
    rcut = np.zeros(npar)
    atoms = centres
    
    change_cutoff_accept = True
    
    while change_cutoff_accept:
        
        print(' Closest approaches           : ')
        for i in range(ntypes):
            for j in range(i,ntypes):
                ic = int(i*(2*ntypes-(i+1))/2+j)
                rcut[ic] = float(input(' pair distance {}-{} : '.format(i+1, j+1)))
        
        rmax = np.max(rcut)  
            
        recalculate_accept = True
        
        while recalculate_accept:
    
            grid = Grid(atoms, vectors)
            grid.make_grid()
            
            nbad = np.zeros(ntypes, dtype=int)
            
            itype = 0
            for i in range(n):
                if i > ncum[itype]:
                    itype += 1
        
                # Determine neighbours 
                
                d, inei = grid.neighbours(i, rmax, 0, n-1)
                
                for j, inj in enumerate(inei):
                    jtype = 0
                    while inj >= ncum[jtype]:
                        jtype += 1
                    i1 = min(itype, jtype)
                    j1 = max(itype, jtype)
                    ic = int(i1*(2*ntypes-(i1+1))/2+j1)
                    if d[j] < rcut[ic]:
                        ind[nbad[itype], itype] = i                
                        nbad[itype] += 1
            
            #print(ind[:nbad[itype],itype])
            #sys.exit()
            totbad = 0
            for i in range(ntypes):
                print('\n {:4d} atoms of type {:2d} have too close neighbours\n'.format(nbad[i], i+1))
                totbad += nbad[i]
        
            # move atom flag
            totbad = 0
            frozen = [False for i in range(ntypes)]
            for i in range(ntypes):
                if nbad[i] == 0:
                    continue
                movatom = input(' Move atoms of type {:d} ? (T/F) : '.format(i+1))
                if strtobool(movatom) == 0:
                    movatom = False
                else:
                    movatom = True
    
                if movatom:
                    frozen[i] = False
                    totbad = totbad + nbad[i]         
            
            move_repeat = True
            if totbad == 0:
                move_repeat = False
                
            # Maximum move
            while move_repeat:
                
                delta = float(input(' Maximum move                 : '))
                #delta = 1.0
                delta = min(delta, vectors[0,0], vectors[1,1], vectors[2,2])
                deltax = delta / np.sqrt(vectors[0,0]**2+vectors[1,0]**2+vectors[2,0]**2)
                deltay = delta / np.sqrt(vectors[0,1]**2+vectors[1,1]**2+vectors[2,1]**2)
                deltaz = delta / np.sqrt(vectors[0,2]**2+vectors[1,2]**2+vectors[2,2]**2)
                niter = int(input(' Max. no. of iterations       : '))
                print()
                    
                #niter = 50000
                ntrys = 0
                
                niter_over = False
                while True:        
                    itype = int(np.random.rand()*ntypes)
                    if frozen[itype] or nbad[itype] == 0: 
                        continue
                    ibad = int(np.random.rand()*nbad[itype])
                    imove = ind[ibad, itype]
                    #print(imove)
                    #if imove == 23:
                    #    print('ibad',ibad)
                    #    sys.exit()
                    xold = atoms[0,imove]
                    yold = atoms[1,imove]
                    zold = atoms[2,imove]
                    dx = 2.*(np.random.rand()-0.5)*deltax
                    dy = 2.*(np.random.rand()-0.5)*deltay
                    dz = 2.*(np.random.rand()-0.5)*deltaz
                    atoms[0,imove] = atoms[0,imove]+dx + 3.
                    atoms[1,imove] = atoms[1,imove]+dy + 3.
                    atoms[2,imove] = atoms[2,imove]+dz + 3.
                    atoms[0,imove] = 2.*(atoms[0,imove]/2.-int(atoms[0,imove]/2.))-1.
                    atoms[1,imove] = 2.*(atoms[1,imove]/2.-int(atoms[1,imove]/2.))-1.
                    atoms[2,imove] = 2.*(atoms[2,imove]/2.-int(atoms[2,imove]/2.))-1.
                    if truncated and  abs(atoms[0,imove])+abs(atoms[1,imove])+abs(atoms[2,imove]) > 1.5:
                       atoms[0,imove] = atoms[0,imove]-sign(1.,atoms[0,imove])
                       atoms[1,imove] = atoms[1,imove]-sign(1.,atoms[1,imove])
                       atoms[2,imove] = atoms[2,imove]-sign(1.,atoms[2,imove])
                      
                    grid.update_grid(imove, xold, yold, zold, atoms[0,imove], atoms[1,imove], atoms[2,imove])
                    
                    d, inei = grid.neighbours(imove, rmax, 0, n-1)
                    
                    accept = True
                    for j, inj in enumerate(inei):
                        jtype = 0
                        while inj >= ncum[jtype]:
                            jtype += 1
                        i1 = min(itype, jtype)
                        j1 = max(itype, jtype)
                        ic = int(i1*(2*ntypes-(i1+1))/2+j1)
                        if d[j] < rcut[ic]:
                            ntrys += 1
                            grid.update_grid(imove, 
                                             atoms[0,imove], atoms[1,imove], atoms[2,imove],
                                             xold, yold, zold)
                            atoms[0,imove] = xold
                            atoms[1,imove] = yold
                            atoms[2,imove] = zold
                            if ntrys >= niter:
                                niter_over = True
                                break
                            accept = False
                            break
                    
                    if niter_over:            
                        break        
                    if not accept: 
                        continue
                    for i in range(ibad+1, nbad[itype]):
                       ind[i-1, itype] = ind[i, itype]        
                    
                    totbad -= 1
                    nbad[itype] = nbad[itype]-1
                    print(" {:4d} atoms of type {:2d} have too close neighbours after {:5d} iterations".format(nbad[itype], itype+1, ntrys))
                    
                    if totbad <= 0:
                        break
                
                if totbad > 0:    
                    move_repeat = bool(input(' Continue ? (T/F)             : '))
                else:
                    move_repeat = False
                        
            re_calculate = input('\n Re-calculate neighbours? (T/F) : ')
            if strtobool(re_calculate) == 0:
                recalculate_accept = False
            else:
                recalculate_accept = True
        
        change_cutoff = input('\n Change cut-offs ? (T/F)      : ')
        if strtobool(change_cutoff) == 0:
            change_cutoff_accept = False
        else:
            change_cutoff_accept = True
    
    return centres

if __name__ == '__main__':
    
    path = './rand.cfg'
    
    ni, centres, vectors = rmc_cfg.read_rmc_config(path)
    
    moveout(ni, centres, vectors)
    nmoltypes = len(ni)
    nsites = np.ones(nmoltypes, dtype=int)
    sites = np.zeros((nmoltypes,3), dtype=int)
    nang = 0    
    truncated = False
    
    rmc_cfg.write_rmc_config('./move.cfg', 'Move configuration', 0, 0, 0, 0, 
                             truncated, vectors, ni, nsites, sites, nang, centres)
    