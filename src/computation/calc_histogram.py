import os
import sys
import math
import numpy as np
from .rmc_configuration import RmcConfiguration

def calc_histogram(cfg, dr):

    ntypes = cfg.nmol_types
    npar = int(ntypes*(ntypes+1)/2)
    ncum = [0]*ntypes
    ncum[0] = cfg.ni[0]
    for itype in range(1,ntypes):
        ncum[itype] = ncum[itype-1]+cfg.ni[itype]

    d = cfg.d
    nr = int(d/dr)+1
    
    # initialize histogram[nr][npar]
    histogram = np.zeros((nr, npar), dtype=int)
    metric = cfg.metric

    type2 = 0
    for i in range(cfg.nmol):
        if(i > ncum[type2]-1):
            type2 = type2 + 1
        type1 = 0
        for j in range(i):
            if(j > ncum[type1]-1):
                type1 = type1 + 1
            ic = int(type1*(2*ntypes-type1-1)/2+type2)
            #ic=(type1-1)*(2*ntypes-type1)/2+type2

            x = cfg.atoms.positions[i][0]-cfg.atoms.positions[j][0]+3.0
            y = cfg.atoms.positions[i][1]-cfg.atoms.positions[j][1]+3.0
            z = cfg.atoms.positions[i][2]-cfg.atoms.positions[j][2]+3.0
            x = 2.0*(x/2.0-int(x/2.0))-1.0
            y = 2.0*(y/2.0-int(y/2.0))-1.0
            z = 2.0*(z/2.0-int(z/2.0))-1.0
            if (cfg.truncated == True and math.fabs(x)+math.fabs(y)+math.fabs(z) > 1.5):
                x = x-sign(1.0, x)
                y = y-sign(1.0, y)
                z = z-sign(1.0, z)
            
            d = metric[0][0]*x*x+metric[1][1]*y*y+metric[2][2]*z*z \
              + 2.0*(metric[0][1]*x*y+metric[0][2]*x*z+metric[1][2]*y*z)

            d = math.sqrt(d)

            # python3 
            # round=lambda x:(x*2+1)//2

            ig=int(round(d/dr))
            if (ig < nr):
                histogram[ig][ic] = histogram[ig][ic]+1

    return histogram

def sign(a, b):
    if(b >= 0.0):
        return math.fabs(a)
    else:
        return -math.fabs(a)


