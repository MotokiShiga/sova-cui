# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 10:51:25 2022

@author: H.Morita
"""

import rmc_configuration
import gridding
import numpy as np
import os

MAX_THETA = 10000

def triplets(nth, nba, files, r, out_file, angcos):
    
    # Read first configuration    
    cfg = rmc_configuration.RmcConfiguration(files[0])
    cfg.read()
    title = cfg.title
    ntypes = len(cfg.ni)
    
    rmax = np.zeros((ntypes,ntypes))
    
    n = 0
    for i in range(ntypes):
        for j in range(i,ntypes):
            rmax[i,j] = r[n]
            n += 1
            
    for i in range(ntypes):
        for j in range(i+1, ntypes):
            rmax[j,i] = rmax[i,j]  
            
    rbig = 0.0
    for i in range(ntypes):
        for j in range(i, ntypes):
            rbig = max(rbig, rmax[i,j])

    nth = min(nth, MAX_THETA)
    dcth = 2.0/nth
    
    cth = [dcth/2.0 + i*dcth-1.0 for i in range(nth)]
    ncth = np.zeros((nth,ntypes,ntypes,ntypes))
        
    for file in files:
        
        if os.path.exists(file) == False:
            print('file not found : ', file)
            return
        
        cfg = rmc_configuration.RmcConfiguration(file)
        cfg.read()
        
        ncum = []
        ncum.append(cfg.ni[0])
        for i in range(1, ntypes):
            ncum.append(ncum[i-1]+cfg.ni[i])
        
        # Set up grid and assign atoms to grid points
        grid = gridding.Grid(cfg)
        grid.make_grid()
        
        metric = cfg.metric
        
        i2 = 0
        for i in range(cfg.nmol):
            if i == ncum[i2]: i2 = i2 + 1
                
            """
            Determine neighbours
            i    : center atom
            rbig : max radius
            0    : begin index
            cfg.nmol-1 : end index
            """
            grid.neighbours(i, rbig, 0, cfg.nmol-1)
            neigh = len(grid.inei)
            inei = grid.inei
            coords = grid.coords
            d = grid.d
        
            if nba > 0:
                neigh = min(neigh, nba)

            itype = np.zeros(neigh, dtype=int)
            for j in range(neigh):                
                for k in range(ntypes):
                    if inei[j] > ncum[k]-1: itype[j] = itype[j]+1
                                
            for j in range(neigh):
                for k in range(j+1, neigh):
                    costh = 0.0
                    for ia in range(3):
                        for ib in range(3):
                            costh = costh + metric[ia][ib] * coords[ia][j] * coords[ib][k]
                
                    costh = costh / d[j] / d[k]
                    costh = max(-1, min(costh, 1.0))
                    ith = int((costh + 1.0) / dcth)
                    i1 = max(itype[j], itype[k])
                    i3 = min(itype[j], itype[k])
                    
                    if d[j] <= rmax[max(itype[j], i2)][min(itype[j], i2)] and \
                       d[k] <= rmax[max(itype[k], i2)][min(itype[k], i2)]:
                           
                        ncth[ith,i1,i2,i3] = ncth[ith,i1,i2,i3] + 1
                        

        # output                        
        cosang = False
        if angcos[0] == 'a' or angcos[0] == 'A': cosang = True

        fw = open(out_file, 'w')

        fw.write('Output from TRIPLETS - bond angle correlations\n')
        fw.write('TITLE\n')
        fw.write(title + '\n')
        fw.write('Maximum r values : \n')
        for i in range(ntypes):
            for j in range(i, ntypes):
                fw.write("%17.8F" % rmax[i][j])
        fw.write('\n\n')

        if cosang == True:
            fw.write('Writing angle distribution\n\n')
        else:
            fw.write('Writing cosine distribution\n\n')

        for i1 in range(ntypes):
            for i2 in range(ntypes):
                for i3 in range(i1+1):
                    fw.write('b{0}{1}{2}\n'.format(i1+1, i2+1, i3+1))
                    fw.write('PLOTS\n')
                    fw.write('%17d' % nth)
                    fw.write('%17d\n' % 1)
                    
                    _sum = 0.
                    for j in range(nth):
                        _sum = _sum + ncth[j,i1,i2,i3]

                    _sum = max(_sum, 1.)
                    for j in range(nth):
                        p = ncth[j,i1,i2,i3] / _sum / dcth
                        if cosang == True:
                            thet = 180. * np.arccos(cth[j])/np.pi
                            fw.write("%17.8F" % thet)                            
                            fw.write("%17.8F\n" % (p*np.sqrt(1.0-cth[j]*cth[j])))
                        else:
                            fw.write("%17.8F" % cth[j])
                            fw.write("%17.8F\n" % p)

                    fw.write('\n')
                                

if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculates triplet correlations')
    parser.add_argument('--nth','-nth', help='No. of theta points', type=int)    
    parser.add_argument('--nba','-nba', help='No. of neighbours for bond ang (0 for all)', default=0)
    parser.add_argument('--files','-files', help='calculate RMC cfg files list',nargs='+')    
    parser.add_argument('--rmax','-rmax', help='Maximum r values (1 1) (1 2) (2 2)...', nargs='+')
    parser.add_argument('--out','-out', help='Output file')
    parser.add_argument('--angcos','-angcos', help='(A)ngle or (C)osine distribution')
    
    args = parser.parse_args()
    
    rmax = [float(n) for n in args.rmax]
        
    triplets(args.nth,
             args.nba, 
             args.files, 
             rmax, 
             args.out, 
             args.angcos)
    
    '''
    # No. of theta pts
    nth = 10
    # No. of neighbours for bond ang (0 for all)
    nba = 0    
    # cfg files   
    files = ['sio.cfg']
        
    out_file = 'sio.tri'
    angcos = 'A'  # A:angle, C:cosine
    '''
    
