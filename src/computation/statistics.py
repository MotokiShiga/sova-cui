#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 12:42:00 2024

@author: morita
"""

import math
import numpy as np
from core.elements import numbers
from core.gridding import metric, d, volume
from . import element_data as elem
import collections

try:
    from computation.histogram import histogram as hist
    has_hist = True
except ImportError:
    import computation.calc_histogram as hist
    has_hist = False
    print('Warning calculate histogram is slow')

def gr(histogram,vectors,numbers,dr):
    _gr = np.zeros_like(histogram, dtype=float)
    nr, npar = histogram.shape  
    _volume = volume(vectors)
    ntypes = len(numbers)
    nxn = np.zeros(npar)
    _r = np.zeros(nr)
    
    ic = 0
    for itype in range(ntypes):
        for jtype in range(itype, ntypes):
            nxn[ic] = numbers[itype]*numbers[jtype]
            #ic = int(itype*(2*ntypes-itype-1)/2+jtype)
            if(itype != jtype):
                nxn[ic] = nxn[ic]*2
            ic += 1
    
    for ir in range(1, nr):
        _r[ir] = ir*dr
        gnorm = (3.0*ir*ir+0.25)*dr*dr*dr*2.0*math.pi/(3.0*_volume)
        for ic in range(npar):
            _gr[ir,ic] = histogram[ir][ic]/(gnorm*nxn[ic])
    
    return _r, _gr

def total_gr(gr,coeff):
    nr, npar = gr.shape
    _total_gr = np.zeros(nr)
    for ir in range(nr):
        for ic in range(npar):
            _total_gr[ir] += coeff[ic]*(gr[ir,ic]-1.0)
    return _total_gr

def histogram(centres,elements,vectors,dr,symbols=None,truncated=False):
    
    def sign(a, b):
        if(b >= 0.0):
            return math.fabs(a)
        else:
            return -math.fabs(a)

    if symbols is None:
        atomic_numbers = numbers
    else:
        atomic_numbers = {str.upper(symbols[i]):(i+1) for i in range(len(symbols))}
    indexes = list(range(len(elements)))
    sorted_elements, sorted_indexes = zip(*sorted(zip(elements, indexes), 
                                                  key=lambda x: (atomic_numbers.get(str.upper(x[0]), x[1])),
                                                  reverse=False))
    counter = collections.Counter(elements)        
    types, ni = zip(*sorted(counter.items(), key=lambda x: atomic_numbers.get(str.upper(x[0]))))
    _d = d(vectors)
    ntypes = len(ni)
    npar = int(ntypes*(ntypes+1)/2)
    types = { types[i] : i for i in range(ntypes)}
    
    if has_hist:
        atoms = []
        for i in sorted_indexes:
            atoms.append(list(centres[i]))
        _metric = []
        for m in metric(vectors):
            _metric.append(list(m))
        ni = list(ni)
        histogram = np.array(hist.calc_histogram(atoms,_metric,ni,_d,dr,truncated))
    else:
        _metric = metric(vectors)
        nr = int(_d/dr)+1
        # initialize histogram[nr][npar]
        histogram = np.zeros((nr, npar), dtype=int)
        
        for i in range(len(elements)):
            type2 = types[elements[i]] 
            for j in range(i):
                type1 = types[elements[j]]
                ic = int(type1*(2*ntypes-type1-1)/2+type2)
                #ic=(type1-1)*(2*ntypes-type1)/2+type2

                x = centres[i][0]-centres[j][0]+3.0
                y = centres[i][1]-centres[j][1]+3.0
                z = centres[i][2]-centres[j][2]+3.0
                x = 2.0*(x/2.0-int(x/2.0))-1.0
                y = 2.0*(y/2.0-int(y/2.0))-1.0
                z = 2.0*(z/2.0-int(z/2.0))-1.0
                if (truncated == True and math.fabs(x)+math.fabs(y)+math.fabs(z) > 1.5):
                    x = x-sign(1.0, x)
                    y = y-sign(1.0, y)
                    z = z-sign(1.0, z)
                
                dis = _metric[0][0]*x*x+_metric[1][1]*y*y+_metric[2][2]*z*z \
                    + 2.0*(_metric[0][1]*x*y+_metric[0][2]*x*z+_metric[1][2]*y*z)

                dis = math.sqrt(dis)

                # python3 
                # round=lambda x:(x*2+1)//2

                ig = int(round(dis/dr))
                if (ig < nr):
                    histogram[ig][ic] = histogram[ig][ic]+1

    return histogram


def ncoeff(symbols, frac, norm=True):
    
    FOUR_PI = 12.5663706143592  # 4*PI

    ntypes = len(symbols)
    npar = int(ntypes*(ntypes+1)/2)
    
    neutron = []
    for s in symbols:
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
            n += 1

    for i in range(npar):
        norm_coeff[i] = coeff[i]/ispars*FOUR_PI
        
    if norm == False:
        return coeff
    else:
        return norm_coeff


def xcoeff(symbols, frac, q, option):
    """
    Options for ceofficients are: 
        1 - Use as calculated
        2 - Divide by <f**2)
        3 - Divide by <f>**2
    """
    
    nq = len(q)
    ntypes = len(symbols)
    npar = int(ntypes*(ntypes+1)/2)
    ff = [0.0]*npar
    xcoeff = [0.0]*nq
    
    for iq in range(nq):    
        q2 = q[iq] * q[iq] / (16.0 * math.pi * math.pi)
        fsq = 0.0
        f = 0.0
        local_xcoeff = [0.0]*npar
        
        for itype in range(ntypes):
            e = symbols[itype]
            a = elem.Xray[e].a
            b = elem.Xray[e].b
            c = elem.Xray[e].c
            
            _ff = 0.0
            for n in range(5):
                _ff = _ff + a[n]*math.exp(-b[n]*q2)
            _ff = _ff + c
            ff[itype] = _ff
            
            f = f + frac[itype]*_ff
            fsq = fsq + frac[itype]*_ff*_ff
         
        for itype in range(ntypes):
            for jtype in range(itype, ntypes):
                ic = int(itype*(2*ntypes-itype-1)/2) + jtype

                local_xcoeff[ic] = frac[itype]*frac[jtype]*ff[itype]*ff[jtype]
                if (itype != jtype):
                    local_xcoeff[ic] = local_xcoeff[ic]*2.0
                if (option == 2):
                    local_xcoeff[ic] = local_xcoeff[ic]/fsq
                if (option == 3):
                    local_xcoeff[ic] = local_xcoeff[ic]/(f*f)
                  
        xcoeff[iq] = local_xcoeff 
        
    return xcoeff
