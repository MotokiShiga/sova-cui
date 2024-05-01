# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:22:04 2023

@author: H. Morita
"""

import numpy as np
from rmc_cfg import write_rmc_config

def rmc_init():
    import datetime
    
    dt = datetime.datetime.now()
    iseed = dt.hour*60*60 + dt.minute*60 + dt.second + dt.microsecond/1000
    iseed = 2*int(iseed)+1    
    return iseed
    
def random(rho, ni, nang=0):
    
    truncated = False
    vectors = np.identity(3)
    nmol = np.sum(ni)
    
    iseed = rmc_init()
    np.random.seed(iseed)
    centres = np.zeros((3+nang,nmol))
    #centres = [[ 0. for i in range(1000)] for j in range(3)]
        
    i = 0
    while i < nmol:
        centres[0,i] = 2.*np.random.rand()-1.
        centres[1,i] = 2.*np.random.rand()-1.
        centres[2,i] = 2.*np.random.rand()-1.
        if nang >= 1: centres[3,i] = np.acos(2.*np.random.rand()-1.)
        if nang >= 2: centres[4,i] = 2.*np.pi*np.random.rand()
        if nang >= 3: centres[5,i] = 2.*np.pi*np.random.rand()
        if truncated and (np.abs(centres[0,i])
                        + np.abs(centres[1,i])
                        + np.abs(centres[2,i]) > 1.5): 
            i -= 1
        i += 1
        
    vectors[0,0] = 0.5*((nmol)/rho)**(1./3.)
    vectors[1,1] = vectors[0,0]
    vectors[2,2] = vectors[0,0]
        
    return vectors, centres
    
if __name__ == '__main__':
    ni = [1000]
    nmoltypes = len(ni)
    vectors, centres = random(0.033, ni)
    nsites = np.ones(nmoltypes, dtype=int)
    sites = np.zeros((nmoltypes,3), dtype=int)
    nang = 0
    
    truncated = False
    
    write_rmc_config('./rand.cfg', 'Random configuration', 0, 0, 0, 0, 
                     truncated, vectors, ni, nsites, sites, nang, centres)
    