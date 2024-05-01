# -*- coding: utf-8 -*-

import numpy as np
import math
from . import element_data as elem
import sys

FOUR_PI = 12.5663706143592  # 4*PI

def calc_ncoeff(symbol, frac, norm=False):

    ntypes = len(symbol)
    npar = int(ntypes*(ntypes+1)/2)
    
    neutron = []
    for s in symbol:
        n = elem.Neutron.first(s)
        neutron.append(n)

    coeff = np.zeros(npar)
    norm_coeff = np.zeros(npar)
    
    n=0
    ispars = 0.0
    for i in range(ntypes):
        for j in range(i,ntypes):
            coeff[n] = frac[i]*frac[j]* \
                       (neutron[i].bc*neutron[j].bc + \
                        neutron[i].c*neutron[j].c) / 100.0
            if i != j:
                coeff[n] = coeff[n]*2.0
                
            ispars = ispars + FOUR_PI*coeff[n]
            n=n+1

    for i in range(npar):
        norm_coeff[i] = coeff[i]/ispars*FOUR_PI
        
    if norm == False:
        return coeff
    else:
        return norm_coeff


if __name__ == '__main__':
    
    symbol = ['Si', 'O']
    frac = [0.33, 0.67]
    norm = True
    
    calc_ncoeff(symbol, frac)
