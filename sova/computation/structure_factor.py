#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 12:42:00 2024

@author: morita
"""

import math
import numpy as np
from ..core.elements import numbers
from ..core.gridding import metric, d, volume
from . import element_data as elem
import collections
from .geometry import Polyhedron

try:
    from .histogram import histogram as hist
    has_hist = True
except ImportError:
    has_hist = False
    print('Warning calculate histogram is slow')

def gr(atoms,histogram,dr):
    numbers = atoms.ni
    vectors = atoms.volume.vectors
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
        gnorm = (3.0*ir*ir+0.25)*dr*dr*dr*2.0*np.pi/(3.0*_volume)
        for ic in range(npar):
            _gr[ir,ic] = histogram[ir][ic]/(gnorm*nxn[ic])
    
    return _r, _gr

def total_gr(gr,coeff):
    nr, npar = gr.shape
    _total_gr = np.zeros(nr)
    for ir in range(nr):
        for ic in range(npar):
            _total_gr[ir] += coeff[ic]*(gr[ir][ic]-1.0)
    return _total_gr + 1.0

def SQ(atoms,gr,qmin,qmax,dr,dq):
    n = atoms.number
    _volume = volume(atoms.volume.vectors)
    nq = int(math.ceil((qmax-qmin)/dq))+1
    q =np.array([(qmin+float(i)*dq) for i in range(nq)])
    nr, npar = gr.shape
    sqr = np.zeros((nr, nq+1), dtype=float)
    _sq = np.zeros((nq, npar), dtype=float)
    
    for iq in range(nq):
        s = np.zeros(npar)
        for ir in range(1, nr):
            r = float(ir)*dr
            sqr = 4.0*np.pi*float(n)/_volume*r*np.sin(r*q[iq])/q[iq]*dr
            for ic in range(npar):
                s[ic] += (gr[ir][ic]-1.0)*sqr
        
        _sq[iq] = s
    
    return q, _sq+1.0

def total_SQ(sq,coeff):
    nq, npar = sq.shape
    total_sq = np.zeros(nq)
    for iq in range(nq):
        for ic in range(npar):
            total_sq[iq] += coeff[ic]*sq[iq,ic]
            
    return total_sq

def total_FQ(sq,coeff):
    nq, npar = sq.shape
    total_fq = np.zeros(nq)
    for i in range(nq):
        for j in range(npar):
            total_fq[i] += (sq[i][j]*coeff[i][j])
    return total_fq

def Gr(r,total_gr,rho):
    Gr = 4.*np.pi*r*rho*total_gr
    return Gr
  
def Tr(r,total_gr,rho):  
    Tr = 4.*np.pi*r*rho*(total_gr+1.0)
    return Tr

def Nr(r,Tr):
    Nr = r*Tr
    return Nr

def histogram(atoms,dr,symbols=None,truncated=False):
    elements = atoms.elements 
    def sign(a, b):
        if(b >= 0.0):
            return math.fabs(a)
        else:
            return -math.fabs(a)
    if symbols is None:
        atomic_numbers = numbers
    else:
        atomic_numbers = {str.upper(symbols[i]):(i+1) for i in range(len(symbols))}
        atoms.symbol_order(symbols)
    indexes = list(range(len(elements)))
    sorted_elements, sorted_indexes = zip(*sorted(zip(elements, indexes), 
                                                  key=lambda x: (atomic_numbers.get(str.upper(x[0]), x[1])),
                                                  reverse=False))
    types, ni = zip(*sorted(atoms.numbers.items(), key=lambda x: atomic_numbers.get(str.upper(x[0]))))
    ntypes = len(ni)
    npar = int(ntypes*(ntypes+1)/2)
    types = { types[i] : i for i in range(ntypes)}
    
    if atoms.volume.periodic:
        centres = atoms.norm_positions
        elements = atoms.elements 
        vectors = atoms.volume.vectors
        _d = d(vectors)
        nr = int(_d/dr)+1
                
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
    else:
        centres = atoms.positions
        elements = atoms.elements 
        length = atoms.maximum - atoms.minimum
        nr = int(min(length)/2/dr)+1
        r = np.zeros(nr)
        histogram = np.zeros((nr, npar), dtype=int)
        for i in range(len(elements)):
            type2 = types[elements[i]] 
            for j in range(i):
                type1 = types[elements[j]]
                ic = int(type1*(2*ntypes-type1-1)/2+type2)
                #ic=(type1-1)*(2*ntypes-type1)/2+type2

                dis = np.linalg.norm(centres[i]-centres[j])

                # python3 
                # round=lambda x:(x*2+1)//2

                ig = int(round(dis/dr))
                if (ig < nr):
                    histogram[ig][ic] = histogram[ig][ic]+1
    
    r = np.zeros(nr)
    for ir in range(1, nr):
        r[ir] = ir*dr

    return r, histogram


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
    
    n = 0
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

def xcoeff(symbols, frac, q, option=3):
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
        q2 = q[iq]*q[iq]/(16.0*np.pi*np.pi)
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

def neighbor(atoms,center,bond_lengths):
    """
    Calculate neighbour indeces and dixtance around center atom.

    Parameters
    ----------
    center : int
        atom index.
    bond_lengths : dict
        maximum bond lengths for atomic element pairs
        Set -1 to pairs for which you don't want to build bonds.
        

    Returns
    -------
    indices : list
        neighbour indices.
    dis : list
        neighbour distances.

    """
    if atoms is None:
        print('Not found Atoms data')
        return
    
    # calculate rmax
    bond_matrix = np.identity(len(atoms.elements_kind))
    for i, elem1 in enumerate(atoms.elements_kind):
        for j, elem2 in enumerate(atoms.elements_kind):
            pair = (elem1, elem2)
            if pair in bond_lengths.keys():
                bond_matrix[i][j] = bond_lengths[pair]
            else:
                pair = (elem2, elem1)
                bond_matrix[i][j] = bond_lengths[pair]
    rmax = np.max(bond_matrix)        
    elements_indices = [atoms.elements_kind.index(e) for e in atoms.elements]
    
    grid = atoms.grid    
    if grid is None:
        position = atoms.positions[center]        
        distance_vector_array = atoms.positions - position
        distance_squared_array = np.sum(np.square(distance_vector_array), axis=1)
        delta_squared = np.square([bond_matrix[elements_indices[center]][j] for j in elements_indices])
        bond_target_indices = (distance_squared_array <= delta_squared).nonzero()[0]
        bond_target_indices = bond_target_indices[bond_target_indices != center]
        inei = bond_target_indices
        dis = np.sqrt(distance_squared_array[bond_target_indices])
        if len(inei):            
            dis, inei =zip(*sorted(zip(dis, inei)))
        inei = np.array(inei)
        dis = np.array(dis)
    else:        
        grid.neighbours(center, rmax, 0, atoms.number)           
        delta = np.array([bond_matrix[elements_indices[center]][elements_indices[j]] for j in grid.inei])
        indices = (grid.d <= delta).nonzero()[0]
        inei = np.array(grid.inei)[indices]
        dis = np.array(grid.d)[indices]
    return inei, dis

def neighbors(atoms,rmin=None,rmax=None):
    """
    Calculate coordination number.

    Parameters
    ----------
    rmin : list, optional
        minimum radius list. pair list order (1,1), (1,2), (2,2).... The default is None.
    rmax : list, optional
        maximum radius list. pair list order (1,1), (1,2), (2,2).... The default is None.
    
    example.
    
    rmin = [0.0, 0.0, 0.0]    
    rmax = [4.0, 2.0, 3.0]    

    Returns
    -------
    None.

    """
    MAX_NEIGH = 20

    if atoms is None:
        print('Not found Atoms data')
        return
    
    grid = atoms.grid
    _rmax = np.max(rmax)
    ntypes = len(atoms.symbols)    
    _neigh = np.zeros((ntypes,ntypes,MAX_NEIGH), dtype=int)
    
    for i in range(atoms.number):            
        grid.neighbours(i, _rmax, 0, atoms.number)
        
        inei = grid.inei
        #icoords = grid.coords
        distances = grid.d

        itype = atoms.symbols.index(atoms.elements[i])
        jtypes = [atoms.symbols.index(atoms.elements[n]) for n in inei]
        
        n = np.zeros(ntypes, dtype=int)
        for jtype, dis in zip(jtypes, distances):
            #if jtype[j] >= itype:
            ic = int(itype*(2*ntypes-(itype+1))/2+jtype)
            if dis <= rmax[ic] and dis >= rmin[ic]:
                n[jtype] += 1
            
        for jtype in range(ntypes):
            _neigh[itype][jtype][n[jtype]-1] += 1
    
    #print('Calculation of neighbours in {}'.format(result.filepath))
    print('Calculation of neighbours')
    print('')
    print('No. of atom types = {}'.format(ntypes))
    print('')
    print('Minimum bond lengths = ', " ".join(list(map(str, rmin))))
    print('Maximum bond lengths = ', " ".join(list(map(str, rmax))))
    print('')
    
    for i in range(ntypes):
        for j in range(ntypes):
            print('{0} - {1} neighbours :\n'.format(atoms.symbols[i],atoms.symbols[j]))
            cdno = 0.
            for k in range(MAX_NEIGH):
               tabline = False
               if _neigh[i,j,k] != 0:
                   tabline = True
                   quot = _neigh[i,j,k]*100./atoms.ni[i]
                   cdno = cdno+quot*(k+1)/100.
               
               if tabline == True:
                   print(' {:>10d} {:>10d} {:>10.3f}'.format(k+1, _neigh[i,j,k], quot))
            
            print('')
            print(' Average coordination {:>6.2f}'.format(cdno))
            print('---------------------------\n')

def coords(atoms, metric, i, j):
    x = atoms.norm_positions[j][0]-atoms.norm_positions[i][0]+3.
    y = atoms.norm_positions[j][1]-atoms.norm_positions[i][1]+3.
    z = atoms.norm_positions[j][2]-atoms.norm_positions[i][2]+3.
    x = 2.*(x/2.-int(x/2.))-1.
    y = 2.*(y/2.-int(y/2.))-1.
    z = 2.*(z/2.-int(z/2.))-1.
    
    d = metric[0][0]*x*x+metric[1][1]*y*y+metric[2][2]*z*z \
      + 2.0*(metric[0][1]*x*y+metric[0][2]*x*z+metric[1][2]*y*z)
    
    return np.array([x,y,z]), math.sqrt(d)

def cosine(a, b, c):
    """
    

    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    b : list 
        center position list.
    c : TYPE
        DESCRIPTION.

    Returns
    -------
    cos_theta : TYPE
        DESCRIPTION.

    """
    ab = np.array(a) - np.array(b)
    bc = np.array(c) - np.array(b)
        
    dot_product = np.dot(ab, bc)
    
    norm_ab = np.linalg.norm(ab)
    norm_bc = np.linalg.norm(bc)
        
    cos_theta = dot_product / norm_ab / norm_bc
    
    return cos_theta

def triplets(atoms,bond_lengths,nth=10,norm_sin=True):
    MAX_THETA = 10000
    
    ntypes = len(atoms.symbols)
    nth = min(nth, MAX_THETA)
    ncth = np.zeros((nth+1,ntypes,ntypes,ntypes), dtype=np.int32)    
    dcth = 180.0/nth
    dth = 180.0/(nth-1)    
    cth = [i*dcth for i in range(nth+1)]

    # calculate angle distribution
    if atoms.grid is not None:
        _metric = metric(atoms.volume.vectors)
        
        for i in range(atoms.number):            
            i2 = atoms.symbols.index(atoms.elements[i])
            inei, dis = neighbor(atoms,i,bond_lengths)
            itypes = [atoms.symbols.index(atoms.elements[n]) for n in inei]
            
            for nj, j in enumerate(inei):
                for nk, k in enumerate(inei):
                    if nj >= nk: continue
                    
                    coords_j, dij = coords(atoms, _metric, i, j)
                    coords_k, dik = coords(atoms, _metric, i, k)                    
                    costh = 0.0
                    for ia in range(3):
                        for ib in range(3):
                            costh = costh + _metric[ia][ib] * coords_j[ia] * coords_k[ib]                    
                    costh = costh / dij / dik
                    costh = max(-1, min(costh, 1.0))
                    thet = 180. * np.arccos(costh)/np.pi
                    ith = int(ceil(thet-dth/2, dth)/dth)
                    i1 = max(itypes[nj], itypes[nk])
                    i3 = min(itypes[nj], itypes[nk])
                    ncth[ith,i1,i2,i3] = ncth[ith,i1,i2,i3] + 1
    else:        
        for i in range(atoms.number):
            i2 = atoms.symbols.index(atoms.elements[i])
            inei, dis = neighbor(atoms,i,bond_lengths)
            itypes = [atoms.symbols.index(atoms.elements[n]) for n in inei]
            
            for nj, j in enumerate(inei):
                for nk, k in enumerate(inei):
                    if nj >= nk: continue

                    costh = cosine(atoms.positions[j], atoms.positions[i], atoms.positions[k])
                    costh = max(-1, min(costh, 1.0))
                    thet = 180. * np.arccos(costh)/np.pi
                    ith = int(ceil(thet-dth/2, dth)/dth)
                    i1 = max(itypes[nj], itypes[nk])
                    i3 = min(itypes[nj], itypes[nk])
                    ncth[ith,i1,i2,i3] = ncth[ith,i1,i2,i3] + 1
        
    _triplets = []
    for i1 in range(ntypes):
        for i2 in range(ntypes):
            for i3 in range(i1+1):                
                _sum = np.sum(ncth[:,i1,i2,i3])
                _triplet = []
                for th, p in zip(cth, ncth[:,i1,i2,i3]):
                    y = p / _sum if _sum != 0.0 else 0.0
                    sin = np.sin(th/180.*np.pi) + 1.0e-6
                    if norm_sin:
                        y /= sin
                    _triplet.append(y)
                _triplets.append(_triplet)
    
    return cth, _triplets
    
    """
    print('Output from TRIPLETS - bond angle correlations')
    print()
    print('Maximum r values : ')
    s = ''
    for i in range(ntypes):
        for j in range(i, ntypes):
            s += "%17.8F" % rmax[i][j]
    print(s)
    print('\n')

    kind = atoms.symbols
    for i1 in range(ntypes):
        for i2 in range(ntypes):
            for i3 in range(i1+1):
                print('{0}-{1}-{2}\n'.format(kind[i1], kind[i2], kind[i3]))
                print('PLOTS')
                print('%17d' % nth + '%17d\n' % 1)                
                
                _sum = np.sum(ncth[:,i1,i2,i3])
                
                for th, p in zip(cth, ncth[:,i1,i2,i3]):
                    y = p / _sum if _sum != 0.0 else 0.0
                    sin = np.sin(th/180.*np.pi) + 1.0e-6
                    if norm_sin == True:
                        y /= sin                
                    print("%17.8F" % th + "%17.8F" % y)
                    
                print()
    """
            
def ceil(x, s):
    return s * math.ceil(float(x)/s)

def floor(x, s):
    return s * math.floor(float(x)/s)

def polyhedra(atoms,center,around,rmax):
    ci = [i for i, s in enumerate(atoms.elements) if s == center]
    ai = [i for i, s in enumerate(atoms.elements) if s == around]
    
    _polyhedra = []
    if atoms.grid is not None:
        grid = atoms.grid
        for i in ci:            
            grid.neighbours(i, rmax, 0, atoms.number-1)
            _neighbors = []
            indexes = []
            for n, nei in enumerate(grid.inei):
                if nei in ai:
                    _neighbors.append(nei)
                    indexes.append(n)            
            coords = np.array(grid.coords).T[indexes]
            vectors = []
            for coord in coords:
                vectors.append(np.dot(coord,atoms.volume.vectors))
            shifts = np.array(grid.shifts).T[indexes]
            p = Polyhedron(atoms,center=i,neighbors=_neighbors,
                           vectors=vectors, shifts=shifts)                        
            _polyhedra.append(p)            
    else:
        # TODO
        pass
    
    return _polyhedra
