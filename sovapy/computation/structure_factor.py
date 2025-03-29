import math
import numpy as np
from ..core import elements
from ..core.gridding import metric, d, volume
from . import element_data as elem
from .geometry import Polyhedron

try:
    from .histogram import histogram as hist
    has_hist = True
except ImportError:
    has_hist = False
    print('Warning calculate histogram is slow')

# def histogram(atoms, dr, symbols=None, truncated=False):
def atom_pair_hist(atoms, dr, truncated=False):
    """Histogram of atom-pair distances
    
    The number of atoms at a distance from r to r+dr from center atoms.

    Parameters
    ----------
    atoms : Atoms object
        Structural data including atomic symbol, position, lattice info,...
    dr : float
        Distance interval (bin width) in the histogram
    trucated : Bool, optional
        Truncated for large distance(?)

    Returns
    -------
    r : numpy.ndarray (#bins)
        List of distances 
    histogram : numpy.ndarray (#r, #atom pairs)
        The number of atoms at distance from r to r+dr from center atoms

    """

    def sign(a, b):
        if(b >= 0.0):
            return math.fabs(a)
        else:
            return -math.fabs(a)
        
    ntypes = len(atoms.num_atoms)
    npar = len(atoms.pairs)
        
    # Dictionary from atomic symbol to the atomic number
    atomic_numbers = elements.numbers
    # List of index
    indexes = list(range(len(atoms.symbols)))
    #Sort atoms (indexes) based on atomic numbers
    sorted_elements, sorted_indexes = zip(*sorted(zip(atoms.symbols, indexes), 
                                                  key=lambda x: (atomic_numbers.get(str.upper(x[0]), x[1])),
                                                  reverse=False))
    # types, ni = zip(*sorted(atoms.num_atoms_dict.items(), key=lambda x: atomic_numbers.get(str.upper(x[0]))))
    #types: atoms.symbol_set, ni: atoms.num_atoms

    # Dictionary from atomic symbol to element index
    symbol2id_dict = { atoms.symbol_set[i] : i for i in range(ntypes)}
    
    # When the structure is periodic
    if atoms.volume.periodic:
        centres = atoms.norm_positions
        vectors = atoms.volume.vectors
        _d = d(vectors)
        nr = int(_d/dr)+1
                
        if has_hist: # Implement using the compiled code
            # Sort atomic positions by sorted indexes
            atom_positions = []
            for i in sorted_indexes:
                atom_positions.append(list(centres[i]))
            
            # Initialize the metric by the lattice vectors
            _metric = []
            for m in metric(vectors):
                _metric.append(list(m))

            # Compute histogram using method compiled from c code.
            histogram = np.array(hist.calc_histogram(atom_positions, _metric, atoms.num_atoms, _d, dr, truncated))

        else: # Implement using the python code

            # Initialize the metric by the lattice vectors
            _metric = metric(vectors)
            # initialize histogram[nr][npar]
            histogram = np.zeros((nr, npar), dtype=int)
            
            # Main loop for computing histograms
            for i in range(len(atoms.symbols)):
                type2 = symbol2id_dict[atoms.symbols[i]] 
                for j in range(i):
                    type1 = symbol2id_dict[atoms.symbols[j]]

                    # Indexing the atom pair
                    if type1 <= type2:
                        ic = int(type1+type2*(1+type2)/2)
                    else:
                        ic = int(type2+type1*(1+type1)/2)
    
                    # Compute the distance between the atom pair
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
    
                    # Update the histogram using the computed distance
                    ig = int(round(dis/dr))
                    if (ig < nr):
                        histogram[ig][ic] = histogram[ig][ic]+1    

    else: # When the structure is non-periodic

        centres = atoms.positions
        length = atoms.maximum - atoms.minimum
        nr = int(min(length)/2/dr)+1
        # r = np.zeros(nr)

        # Main loop for computing histograms
        histogram = np.zeros((nr, npar), dtype=int)
        for i in range(len(atoms.symbols)):
            type2 = symbol2id_dict[atoms.symbols[i]] 
            for j in range(i):
                type1 = symbol2id_dict[atoms.symbols[j]]
                # Indexing atom pair
                if type1 <= type2:
                    ic = int(type1+type2*(1+type2)/2)
                else:
                    ic = int(type2+type1*(1+type1)/2)

                # Compute the distance between the atom pair
                dis = np.linalg.norm(centres[i]-centres[j])

                # Update the histogram using the computed distance
                ig = int(round(dis/dr))
                if (ig < nr):
                    histogram[ig][ic] = histogram[ig][ic]+1
    
    # Generate the list of distance values
    r = np.zeros(nr)
    for ir in range(1, nr):
        r[ir] = ir*dr

    return r, histogram

# def gr(atoms, histogram, dr):
def pair_dist_func(atoms, r, histogram):
    """Pair distribution function (PDF)

    Parameters
    ----------
    atoms : Atoms object
        Structural data including atomic symbol, position, lattice info,...
    r : numpy.ndarray (#bins)
        List of distances
    histogram : numpy.ndarray (#r, #atom pairs) 
        The number of atoms at distance from r to r+dr from the center atom

    Returns
    -------
    partial_gr : numpy.ndarray (#r, #atom pairs)
        Pair distribution functions
    """

    num_atoms = atoms.num_atoms
    vectors = atoms.volume.vectors
    partial_gr = np.zeros_like(histogram, dtype=float)
    nr, npar = histogram.shape  
    _volume = volume(vectors)
    ntypes = len(num_atoms)
    nxn = np.zeros(npar)
    dr = r[1]-r[0]
    
    # Compute products of the numbers of atoms
    for type1 in range(ntypes):
        for type2 in range(type1, ntypes):
            ic = int(type1+type2*(1+type2)/2)
            nxn[ic] = num_atoms[type1]*num_atoms[type2]
            if(type1 != type2):
                nxn[ic] = nxn[ic]*2
    
    # Compute PDF function
    for ir in range(1, nr):
        gnorm = r[ir]*r[ir]*dr*2.0*np.pi/(_volume)
        for ic in range(npar):
            partial_gr[ir,ic] = histogram[ir][ic]/(gnorm*nxn[ic])
    
    return partial_gr

# def SQ(atoms,gr,qmin,qmax,dr,dq):
def partial_structure_factor(atoms, r, partial_gr, qmin, qmax, dq):
    """Partial structure factors

    Parameters
    ----------
    atoms : Atoms object
        Structural data including atomic symbol, position, lattice info,...
    r : numpy.ndarray (#bins)
        List of distances
    partial_gr : numpy.ndarray (#r, #atom pairs)
        Pair distribution functions
    qmin : float
        Minimum value of q (distance in reciprocal space)
    qmax : float
        Maximum value of q (distance in reciprocal space)
    qr : float
        Reciprocal distance interval (bin width) in the histogram

    Returns
    -------
    q : numpy.ndarray (#bins)
        List of reciprocal distances
    partial_sq : numpy.ndarray (#q, #atom pairs)
        Pair structure factors
    """

    # Coefficient
    rho = atoms.atom_number_density
    const = 4.0*np.pi*rho

    #The number of q values
    nq = int(math.ceil((qmax-qmin)/dq))+1
    # q values
    q =np.array([(qmin+float(i)*dq) for i in range(nq)])
    # The number of r values, the number of atom pairs
    nr, npar = partial_gr.shape 
    # sqr = np.zeros((nr, nq+1), dtype=float)
    partial_sq = np.zeros((nq, npar), dtype=float)
    dr = r[1]-r[0]
    
    for iq in range(nq):
        s = np.zeros(npar)
        for ir in range(1, nr):
            # r = float(ir)*dr
            sqr = np.sin(r[ir]*q[iq])/q[iq]*r[ir]*dr
            for ic in range(npar):
                s[ic] += (partial_gr[ir][ic]-1.0)*sqr
        partial_sq[iq,:] = s
    partial_sq = 1.0 + const*partial_sq
    
    return q, partial_sq

# def total_gr(gr,coeff):
def atomc_pair_dist_func_neutron(partial_gr, coeff_neutron):
    """ Atomic pair distribution function by Neutron
    
    Atomic pair distribution functioncomputed by neutron scattaring length 

    Parameters
    ----------
    partial_gr : numpy.ndarray (#r, #atom pairs)
        Pair distribution functions
    coeff_neutron: numpy.ndarray (#atom pairs)
        Coefficients regarding neutron scattering

    Returns
    -------
    total_gr : numpy.ndarray (#r)
        Atomic pair distribution function by neutron diffraction
    """

    # The number of r values, the number of atom pairs
    nr, npar = partial_gr.shape
    total_gr = np.zeros(nr)
    for ir in range(nr):
        for ic in range(npar):
            total_gr[ir] += coeff_neutron[ic]*(partial_gr[ir][ic]-1.0)
    return total_gr + 1.0

# def total_SQ(sq, coeff):
def structure_factor_neutron(partial_sq, coeff_neutron):
    """Structure factor by neutron diffraction 

    Parameters
    ----------
    partial_sq : numpy.ndarray (#q, #atom pairs)
        pair structure factors
    coeff_neutron: numpy.ndarray (#atom pairs)
        Coefficients regarding neutron scattering

    Returns
    -------
    sq_neutron : numpy.ndarray (#q)
        Structure factor by neutron diffraction 
    """

    num_q, num_pairs = partial_sq.shape
    sq_neutron = np.zeros(num_q)
    for i_q in range(num_q):
        for i_pair in range(num_pairs):
            sq_neutron[i_q] += coeff_neutron[i_pair]*partial_sq[i_q, i_pair]
            
    return sq_neutron

# def total_FQ(sq, coeff):
def structure_factor_xray(partial_sq, coeff_xray):
    """Structure factor by X-ray diffraction 

    Parameters
    ----------
    partial_sq : anumpy.ndarray (#q, #atom pairs)
        Pair structure factors
    coeff_xray: numpy.ndarray (#q, #atom pairs)
        Coefficients regarding X-ray scattering

    Returns
    -------
    sq_xray : numpy.ndarray (#q)
        Structure factor by X-ray diffraction 
    """

    num_q, num_pairs = partial_sq.shape
    sq_xray = np.zeros(num_q)
    for i_q in range(num_q):
        for i_pair in range(num_pairs):
            sq_xray[i_q] += (coeff_xray[i_q, i_pair]*partial_sq[i_q, i_pair])
    return sq_xray

# def Gr(r, total_gr, rho):
def reduced_pair_dist_func(r, total_gr, rho):
    """Reduced atomic pair distribution function 
    
    Parameters
    ----------
    r : numpy.ndarray (#bins)
        List of real space distances
    total_gr : numpy.ndarray (#r)
        Atomic pair distribution function by neutron diffraction
    rho : float
        Atomic number density in Angstrom^{-3}

    Returns
    -------
    Gr : numpy.ndarray (#r)
        Reduced atomic pair distribution function
    """

    Gr = 4.*np.pi*r*rho*(total_gr-1)
    return Gr
  
# def Tr(r, total_gr, rho):
def total_corr_fun(r, total_gr, rho):  
    """Total correlation function 
    
    Parameters
    ----------
    r : numpy.ndarray (#bins)
        List of real space distances
    total_gr : numpy.ndarray (#r)
        Atomic pair distribution function by neutron diffraction
    rho : float
        Atomic number density in Angstrom^{-3}

    Returns
    -------
    Tr : numpy.ndarray (#r)
        Total correlation function
    """

    Tr = 4.*np.pi*r*rho*total_gr
    return Tr

# def Nr(r, Tr):
def radial_dist_fun(r, total_gr, rho):
    """Radial distribution function
    
    Parameters
    ----------
    r : numpy.ndarray (#bins)
        List of real space distances
    total_gr : numpy.ndarray (#r)
        Atomic pair distribution function by neutron diffraction
    rho : float
        Atomic number density in Angstrom^{-3}

    Returns
    -------
    Nr : numpy.ndarray (#r)
        Radial distribution function
    """

    Tr = total_corr_fun(r, total_gr, rho)
    Nr = r*Tr
    return Nr

def ncoeff(atoms):
    """Coefficients regarding Neutron diffraction

    Parameters
    ----------
    atoms : data.Atoms object
        Structural data

    Returns
    -------
    coeff : numpy.ndarray (#atom pairs)
        Coefficients regarding Neutron diffraction
    """
    
    npar = len(atoms.pairs)
    
    # Get data containing neutron scattering lengths
    neutron = dict()
    for s0,s1 in atoms.pairs:
        if s0 not in neutron.keys():
            neutron[s0] = elem.Neutron.first(s0)
        if s1 not in neutron.keys():
            neutron[s1] = elem.Neutron.first(s1)
        
    # Compute coefficients of all atom pairs
    n, norm_term = 0, 0.0
    coeff = np.zeros(npar)
    for s0,s1 in atoms.pairs:
        frac0 = atoms.num_atoms_dict[s0]/atoms.num_total
        frac1 = atoms.num_atoms_dict[s1]/atoms.num_total
        # coeff[n] = frac0*frac1* \
        #             (neutron[s0].bc*neutron[s1].bc + \
        #             neutron[s0].c*neutron[s1].c) / 100.0
        coeff[n] = frac0 * frac1 * neutron[s0].bc * neutron[s1].bc
        if s0 != s1:
            coeff[n] = coeff[n]*2.0
        norm_term = norm_term + coeff[n]
        n += 1

    # Normalization
    for i in range(npar):
        coeff[i] = coeff[i] / norm_term
        
    return coeff

def xcoeff(atoms, q):
    """Coefficients regarding X-ray diffraction

    Parameters
    ----------
    atoms : data.Atoms object
        Structural data
    q : numpy.ndarray (#q)
        List of reciprocal distances
    
    Returns
    -------
    coeff_xray: numpy.ndarray (#q, #atom pairs)
        Coefficients regarding X-ray diffraction
    """
    
    num_q = q.size
    ntypes = len(atoms.symbol_set)
    num_pairs = len(atoms.pairs)

    ff = np.zeros(num_pairs)
    coeff_xray = np.zeros((num_q,num_pairs))

    # dictionary from atomic symbol to element index
    atom_id_dict = { atoms.symbol_set[i] : i for i in range(ntypes)}
    frac = atoms.frac_atoms
    
    for iq in range(num_q):
        q2 = q[iq]*q[iq]/(16.0*np.pi*np.pi)
        f = 0.0
        local_coeff_xray = np.zeros(num_pairs)
        
        for itype in range(ntypes):
            # Get coefficients of scattering factor
            s = atoms.symbol_set[itype]
            a = elem.Xray[s].a
            b = elem.Xray[s].b
            c = elem.Xray[s].c
            
            # Compute scattering factor
            _ff = 0.0
            for n in range(5):
                _ff = _ff + a[n]*math.exp(-b[n]*q2)
            _ff = _ff + c
            ff[itype] = _ff
    
            # Add to the normalization term
            f = f + frac[itype]*_ff
         
        # Compute coefficients for each atom pair
        ic = 0
        for s0,s1 in atoms.pairs:
            i0 = atom_id_dict[s0]
            i1 = atom_id_dict[s1]
            local_coeff_xray[ic] = frac[i0]*frac[i1]*ff[i0]*ff[i1]
            if (s0 != s1):
                local_coeff_xray[ic] = local_coeff_xray[ic]*2.0
            local_coeff_xray[ic] = local_coeff_xray[ic] / (f*f)
            ic = ic + 1
                  
        coeff_xray[iq,:] = local_coeff_xray 
        
    return coeff_xray

def neighbor(atoms, center, bond_lengths):
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
    bond_matrix = np.identity(len(atoms.symbol_set))
    for i, elem1 in enumerate(atoms.symbol_set):
        for j, elem2 in enumerate(atoms.symbol_set):
            pair = (elem1, elem2)
            if pair in bond_lengths.keys():
                bond_matrix[i][j] = bond_lengths[pair]
            else:
                pair = (elem2, elem1)
                bond_matrix[i][j] = bond_lengths[pair]
    rmax = np.max(bond_matrix)        
    elements_indices = [atoms.symbol_set.index(e) for e in atoms.symbols]
    
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
        grid.neighbours(center, rmax, 0, atoms.num_total)           
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
    ntypes = len(atoms.symbol_set)    
    _neigh = np.zeros((ntypes,ntypes,MAX_NEIGH), dtype=int)
    
    for i in range(atoms.num_total):            
        grid.neighbours(i, _rmax, 0, atoms.num_total)
        
        inei = grid.inei
        #icoords = grid.coords
        distances = grid.d

        itype = atoms.symbol_set.index(atoms.symbols[i])
        jtypes = [atoms.symbol_set.index(atoms.symbols[n]) for n in inei]
        
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
            print('{0} - {1} neighbours :\n'.format(atoms.symbol_set[i],atoms.symbol_set[j]))
            cdno = 0.
            for k in range(MAX_NEIGH):
               tabline = False
               if _neigh[i,j,k] != 0:
                   tabline = True
                   quot = _neigh[i,j,k]*100./atoms.num_atoms[i]
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

# def triplets(atoms, num_bins=10, norm_sin=True, prob=True, bond_lengths=None):
def bond_angle_hist(atoms, num_bins=10, norm_sin=True, prob=True, bond_lengths=None):
    """Bond angle histograms

    Coompute bond angle histograms using bond lengths setting

    Parameters
    ----------
    atoms : core.Atoms object
        Structural data
    num_bins : int, optional
        The number of bins, by default 10 (Old name: nth)
    norm_sin : bool, optional
        If the histogram is normalize by sin, by default True
    prob : bool, optional
        If the histogram is normalize by the total counts, by default True
    bond_lengths : Atoms.bond_lengths, optional
        Bond length setting, by default None

    Returns
    -------
    angles : numpy.ndarray
        List of angles (Degree) (Old name: cth)
    hist_angles : List (#trios, #bins)
        List of bond angle histograms (Old name: _triplets)
        The order of atom triplets is identical to Atoms.trios
    """
    MAX_THETA = 10000
    
    if bond_lengths is None:
        bond_lengths = atoms.bond_lengths

    # The number of elements (atoms)
    ntypes = len(atoms.symbol_set)
    # The number of histogram bins
    num_bins = min(num_bins, MAX_THETA)
    # Histogram of angles (temporal)
    ncth = np.zeros((num_bins+1, ntypes, ntypes, ntypes), dtype=np.int32)    
    # Interval of angles (the width of histogram bins) (Old nane: dth)
    delta_angle = 180.0/(num_bins-1)    
    
    # dcth = 180.0/num_bins
    #angles = [i*dcth for i in range(num_bins+1)]
    angles = np.linspace(0, 180, num_bins)

    # calculate angle distribution
    if atoms.grid is not None:
        _metric = metric(atoms.volume.vectors)
        
        for i in range(atoms.num_total):            
            i2 = atoms.symbol_set.index(atoms.symbols[i])
            inei, dis = neighbor(atoms, i, bond_lengths)
            itypes = [atoms.symbol_set.index(atoms.symbols[n]) for n in inei]
            
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
                    ith = int(ceil(thet-delta_angle/2, delta_angle)/delta_angle)
                    i1 = max(itypes[nj], itypes[nk])
                    i3 = min(itypes[nj], itypes[nk])
                    ncth[ith,i1,i2,i3] = ncth[ith,i1,i2,i3] + 1
    else:        
        for i in range(atoms.num_total):
            i2 = atoms.symbol_set.index(atoms.symbols[i])
            inei, dis = neighbor(atoms, i, bond_lengths)
            itypes = [atoms.symbol_set.index(atoms.symbols[n]) for n in inei]
            
            for nj, j in enumerate(inei):
                for nk, k in enumerate(inei):
                    if nj >= nk: continue

                    costh = cosine(atoms.positions[j], atoms.positions[i], atoms.positions[k])
                    costh = max(-1, min(costh, 1.0))
                    thet = 180. * np.arccos(costh)/np.pi
                    ith = int(ceil(thet-delta_angle/2, delta_angle)/delta_angle)
                    i1 = max(itypes[nj], itypes[nk])
                    i3 = min(itypes[nj], itypes[nk])
                    ncth[ith,i1,i2,i3] = ncth[ith,i1,i2,i3] + 1
        
    # Histogram of bond angles for outputs
    hist_angles = []
    for i1 in range(ntypes):
        for i2 in range(ntypes):
            for i3 in range(i1+1):                
                _sum = np.sum(ncth[:,i1,i2,i3])
                hist = []
                for th, p in zip(angles, ncth[:,i1,i2,i3]):
                    if prob: # Convert to probability
                        y = p / _sum if _sum != 0.0 else 0.0
                    else:
                        y = p
                    if norm_sin:
                        y /= (np.sin(th/180.*np.pi) + 1.0e-6)
                    hist.append(y)
                hist_angles.append(hist)
    
    return angles, hist_angles
            
def ceil(x, s):
    return s * math.ceil(float(x)/s)

def floor(x, s):
    return s * math.floor(float(x)/s)

# def polyhedra(atoms, center, around, rmax):
def tetrahedra(atoms, center, around, rmax):
    ci = [i for i, s in enumerate(atoms.symbols) if s == center]
    ai = [i for i, s in enumerate(atoms.symbols) if s == around]
    
    _tetrahedra = []
    if atoms.grid is not None:
        grid = atoms.grid
        for i in ci:            
            grid.neighbours(i, rmax, 0, atoms.num_total-1)
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
            p = Polyhedron(atoms, center=i, neighbors=_neighbors,
                           vectors=vectors, shifts=shifts)                        
            _tetrahedra.append(p)            
    else:
        # TODO
        pass
    
    return _tetrahedra
