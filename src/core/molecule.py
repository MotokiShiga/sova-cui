import numpy as np
from numpy import linalg
from scipy import stats
# from scipy.stats import gaussian_kde
import time, itertools, os, shutil, pickle, multiprocessing
import numba as nb
import networkx as nx
from core.gridding import Grid

try:
    import igraph as ig
    has_igraph = True
except ImportError:
    has_igraph = False

"""
_progress = None
_finish = None

def progress(*args):
    if callable(_progress):
        return _progress(*args)

def finish(*args):
    if callable(_finish):
        return _finish(*args)
"""
    
'''
def set_output_callbacks(progress_func, print_func, finished_func, error_func, log_func):
    global _progress, _print_message, _finish, _error, _log

    _progress = progress_func
    _print_message = print_func
    _finish = finished_func


def set_output_callbacks(progress_func, finished_func):
    global _progress, _finish

    _progress = progress_func
    _finish = finished_func
'''    

@nb.jit('Tuple((float32[:], uint32[:,:]))(float32[:,:], float32, float32)')
def jit_distance(xyz, r_max, p_pair=np.float32(0.3)):
    """calculate distances which are shoter than threshold
    
    Args:
        xyz (float[:,:]): 3 dimension coordinates of atoms
        r_max (float): threshold on distance calculations
        p_pair (float): the rate of the number of atom pairs
    
    Returns:
        distance (float[:]): distances between atom pairs
        index_atoms (int[:,:]): ids of atom pairs
    """

    r_max2 = r_max**2

    num_atoms = xyz.shape[0]
    num_pairs = int(num_atoms*(num_atoms-1)/2)
    if num_atoms>1000:
        num_pairs = int(num_pairs*p_pair)
    distance = np.zeros(num_pairs, dtype=np.float32)
    index_atoms = np.zeros((num_pairs,2), dtype=np.uint32)
    cnt = 0
    for i in range(num_atoms):
        for j in range(i):
            d = np.sum((xyz[i,:] - xyz[j,:])**2)
            if (d<r_max2):
                distance[cnt] = np.sqrt(d)
                index_atoms[cnt,:] = [i, j]
                cnt += 1
    distance = distance[:cnt]
    index_atoms = index_atoms[:cnt, :]
        
    return distance, index_atoms

@nb.jit('Tuple((float32[:], uint32[:,:], int32[:,:]))(float32[:,:], float32, float32, float32)')
def jit_distance_supercell(xyz, l_size, r_max, p_pair=np.float32(0.3)):
    """calculate distances which are shoter than threshold
    
    Args:
        xyz (float[:,:]): 3 dimension coordinates of atoms
        l_size  (float): the size of lattice (cube)
        r_max (float): threshold on distance calculations
        p_pair (float): the rate of the number of atom pairs
    
    Returns:
        distance (float[:]): distances between atom pairs
        index_atoms (int[:,:]): ids of atom pairs
        index_cells (int[:,:]): ids of super cell
    """

    r_max2 = r_max**2

    num_atoms = xyz.shape[0]
    num_pairs = int(num_atoms*(num_atoms-1)/2)
    if num_atoms>1000:
        num_pairs = int(num_pairs*p_pair)
    distance = np.zeros(num_pairs, dtype=np.float32)
    index_atoms = np.zeros((num_pairs,2), dtype=np.uint32)
    index_cells = np.zeros((num_pairs,3), dtype=np.int32)
    cnt = 0
    
    for i in range(num_atoms):
        for j in range(i):
            for sx in [-1,0,1]:
                for sy in [-1,0,1]:
                    for sz in [-1,0,1]:
                        s = np.array([sx,sy,sz])*l_size
                        d = np.sum((xyz[i,:] - (xyz[j,:] + s))**2)
                        if (d<r_max2):
                            distance[cnt] = np.sqrt(d)
                            index_atoms[cnt,:] = [i, j]
                            index_cells[cnt,:] = [sx, sy, sz]
                            cnt += 1
    distance = distance[:cnt]
    index_atoms = index_atoms[:cnt, :]
    index_cells = index_cells[:cnt, :]
        
    return distance, index_atoms, index_cells

@nb.jit('Tuple((float32[:], uint32[:,:], int32[:,:]))(float32[:,:], float32[:,:], float32, float32)')
def distance_supercell_noncube(xyz, lattice_matrix, r_max, p_pair=np.float32(1.0)):
    """calculate distances which are shoter than threshold
    
    Args:
        xyz (float[:,:]): 3 dimension coordinates of atoms
        lattice_matrix  (float[:,:]): matrix to convert fractional coordinate to cartesian coordinate
        r_max (float): threshold on distance calculations
        p_pair (float): the rate of the number of atom pairs
    
    Returns:
        distance (float[:]): distances between atom pairs
        index_atoms (int[:,:]): ids of atom pairs
        index_cells (int[:,:]): ids of super cell
    """

    # convert to contiguous arrays in memory (C order).
    xyz = np.ascontiguousarray(xyz)
    lattice_matrix = np.ascontiguousarray(lattice_matrix)

    r_max2 = r_max**2

    num_atoms = xyz.shape[0]
    num_pairs = int(num_atoms*(num_atoms-1)/2)
    if num_atoms>1000:
        num_pairs = int(num_pairs*p_pair)
    distance = np.zeros(num_pairs, dtype=np.float32)
    index_atoms = np.zeros((num_pairs,2), dtype=np.uint32)
    index_cells = np.zeros((num_pairs,3), dtype=np.int32)
    cnt = 0
    
    xyz_frac = xyz @ np.linalg.inv(lattice_matrix)
    
    for i in range(num_atoms):
        for j in range(i):
            for sx in [-1,0,1]:
                for sy in [-1,0,1]:
                    for sz in [-1,0,1]:
                        s = np.array([sx,sy,sz])
                        d = np.sum((xyz[i,:] - (xyz_frac[j,:] + s)@lattice_matrix)**2)
                        if (d<r_max2):
                            distance[cnt] = np.sqrt(d)
                            index_atoms[cnt,:] = [i, j]
                            index_cells[cnt,:] = [sx, sy, sz]
                            cnt += 1
                        # if cnt==(num_pairs-1):
                        #     distance = np.r_[distance, np.zeros(num_pairs, dtype=np.float32)]
                        #     index_atoms = np.r_[index_atoms, np.zeros((num_pairs,2), dtype=np.uint32)]
                        #     index_cells = np.r_[index_cells, np.zeros((num_pairs,3), dtype=np.int32)]
                        #     num_pairs = num_pairs*2
    distance = distance[:cnt]
    index_atoms = index_atoms[:cnt, :]
    index_cells = index_cells[:cnt, :]
        
    return distance, index_atoms, index_cells

def distance_supercell(xyz, lattice_size, r_max, p_pair=0.3):
    """calculate distances which are shoter than threshold
    
    Args:
        xyz (float[:,:]): 3 dimension coordinates of atoms
        lattice_size  (float): the size of lattice (cube)
        r_max (float): threshold on distance calculations
        p_pair (float): the rate of the number of atom pairs
    
    Returns:
        distance (float[:]): distances between atom pairs
        index_atoms (int[:,:]): ids of atom pairs
        index_cells (int[:,:]): ids of super cell
    """

    lattice_size = np.array(lattice_size)
    if lattice_size.size==1:
        lattice_size = np.ones((1,3))*lattice_size
    elif lattice_size.size==3:
        lattice_size = lattice_size.reshape(1,3)
    else:
        print('The size of lattice_size is wrong!')
        return -1

    r_max2 = r_max**2
    num_atoms = xyz.shape[0]
    num_pairs = int(num_atoms*(num_atoms-1)/2)
    if num_atoms>1000:
        num_pairs = int(num_pairs*p_pair)
    
    d_xyz = r_max # length of sub-cell
        
    xyz = xyz - np.min(xyz, axis=0, keepdims=True)
    xyz = xyz % lattice_size
    
    i_xyz = (xyz/d_xyz)
    i_xyz = i_xyz.astype(np.int)

    N = (lattice_size/d_xyz).astype(np.int)+1
    N = N.flatten()
    
    i_sublattice = list()
    for _ in range(N[0]):
        tmp2 = list()
        for _ in range(N[1]):
            tmp1 = list()
            for _ in range(N[2]):
                tmp1.append(list([]))
            tmp2.append(tmp1)
        i_sublattice.append(tmp2)
    for i in range(num_atoms):
        i_sublattice[i_xyz[i,0]][i_xyz[i,1]][i_xyz[i,2]] += [i]
    
    distance = np.zeros(num_pairs, dtype=np.float32)
    index_atoms = np.zeros((num_pairs,2), dtype=np.uint32)
    index_cells = np.zeros((num_pairs,3), dtype=np.int32)
    
    #progress(0)
    cnt = 0
    for nx1 in range(N[0]):
        for ny1 in range(N[1]):
            for nz1 in range(N[2]):
                for sx in [-2, -1, 0, 1, 2]:
                    for sy in [-2, -1, 0, 1, 2]:
                        for sz in [-2, -1, 0, 1, 2]:
                            nx2 = (nx1+sx)%N[0]
                            ny2 = (ny1+sy)%N[1]
                            nz2 = (nz1+sz)%N[2]

                            # position of cell of atom2
                            c_xyz = np.floor(np.array([(nx1+sx)/N[0],(ny1+sy)/N[1],(nz1+sz)/N[2]])).astype(np.int)

                            set1 = i_sublattice[nx1][ny1][nz1]
                            set2 = i_sublattice[nx2][ny2][nz2]
                            
                            for atom1 in set1:
                                for atom2 in set2:
                                    if atom1<atom2:
                                        r2 = np.sum((xyz[atom1,:]-(xyz[atom2,:]+c_xyz*lattice_size.flatten()))**2)         
                                        if (r2<r_max2):
                                            distance[cnt] = np.sqrt(r2)
                                            index_atoms[cnt,:] = [atom1,atom2]
                                            index_cells[cnt,:] = c_xyz
                                            cnt += 1
        #progress(int((nx1+1)/N[0]*100))

    distance = distance[:cnt]
    index_atoms = index_atoms[:cnt,:]
    index_cells = index_cells[:cnt,:]
        
    return distance, index_atoms, index_cells


def enumerate_primitive_ring(atoms_extracted, atoms_all, chemical_bond_index, cutoff_size):
    """Enumerate primitive rings
    
    Parameters
    ----------
    num_atoms : int
        the number of atoms
    chemical_bond_index : list
        the list of atom pairs with bonded
    cutoff_size : int
        the maximum number of atoms in enumerated rings
    
    Returns
    -------
    set_rings_all: set (list of int)
        the enumerated primitive rings
    """

    # initialize graph object
    G = nx.Graph()
    G.add_nodes_from(atoms_all)
    G.add_edges_from(chemical_bond_index)
    D = dict(nx.all_pairs_shortest_path_length(G, cutoff=round(cutoff_size/2)+1))

    if has_igraph:        
        G_ig = ig.Graph()
        G_ig.add_vertices(len(atoms_all))
        G_ig.add_edges(chemical_bond_index)

    set_rings=set()
    #progress(0)
    num = len(atoms_extracted)
    for ic, n_source in enumerate(atoms_extracted):
        length = D[n_source]
        for L in range(1,round(cutoff_size/2)+1):
            i = np.array(list(length.values()))==L
            node_check = np.array(list(length.keys()))[i]
            
            # rings with even number of nodes
            for n_target in node_check:
                if has_igraph:
                    paths = [p for p in G_ig.get_all_shortest_paths(n_source,n_target)]
                else:
                    paths = [p for p in nx.all_shortest_paths(G,source=n_source,target=n_target)]                
                if (len(paths)>1)&(n_source<n_target):
                    #enumerate pairs of paths
                    for p0, p1 in itertools.combinations(paths,2):
                        if len(set(p0)&set(p1))==2: #common nodes are only n_source and n_target
                            ring_tmp = nx.Graph()
                            nx.add_path(ring_tmp, p0)
                            nx.add_path(ring_tmp, p1[::-1])

                            # flag = True
                            # set_node = p0 + p1[-2:0:-1]
                            # L_mid = round(len(set_node)/2)-1
                            # for i in range(len(set_node)):
                            #     j = (i + L_mid)%len(set_node)
                            #     ni, nj = set_node[i], set_node[j]
                            #     if D[ni][nj] < nx.shortest_path_length(ring_tmp, source=ni, target=nj):
                            #         flag=False
                            #         break

                            flag = True
                            set_node = p0 + p1[-2:0:-1]
                            for ni,nj in  itertools.combinations(set_node,2):
                                if D[ni][nj] < nx.shortest_path_length(ring_tmp, source=ni, target=nj):
                                    flag=False
                                    break

                            if flag:
                                # to arrange order
                                path = p0 + p1[-2:0:-1]
                                path = np.array(path)
                                i_min = np.argmin(path)
                                path = tuple(np.r_[path[i_min:],path[:i_min]])
                                # to aline the orientation
                                if path[-1]<path[1]:
                                    path = path[::-1]
                                    path = tuple(np.r_[path[-1],path[:-1]])
                                if path in set_rings:
                                    pass
                                else:
                                    set_rings.add(path)

            # rings with odd number of nodes
            for nt0, nt1 in itertools.combinations(node_check,2):
                if nt0 in G.neighbors(nt1): # if nt0 and nt1 are linked
                    if has_igraph:
                        paths0 = [p for p in G_ig.get_all_shortest_paths(n_source,nt0)]
                        paths1 = [p for p in G_ig.get_all_shortest_paths(n_source,nt1)]
                    else:
                        paths0 = [p for p in nx.all_shortest_paths(G,source=n_source,target=nt0)]
                        paths1 = [p for p in nx.all_shortest_paths(G,source=n_source,target=nt1)]
                    for p0,p1 in itertools.product(paths0, paths1):
                        if len(set(p0)&set(p1))==2: #common nodes are only n_source and n_target
                            ring_tmp = nx.Graph()
                            nx.add_path(ring_tmp, p0)
                            nx.add_path(ring_tmp, [nt0,nt1])
                            nx.add_path(ring_tmp, p1[::-1])

                            # flag = True
                            # set_node = p0 + p1[-2:0:-1]
                            # L_mid = round(len(set_node)/2)-1
                            # for i in range(len(set_node)):
                            #     j = (i + L_mid)%len(set_node)
                            #     ni, nj = set_node[i], set_node[j]
                            #     if D[ni][nj] < nx.shortest_path_length(ring_tmp, source=ni, target=nj):
                            #         flag=False
                            #         break

                            flag = True
                            set_node = p0 + p1[-2:0:-1]
                            for ni,nj in  itertools.combinations(set_node,2):
                                if D[ni][nj] < nx.shortest_path_length(ring_tmp, source=ni, target=nj):
                                    flag=False
                                    break

                            if flag:
                                # to arrange order
                                path = p0 + p1[-1:0:-1]
                                path = np.array(path)
                                i_min = np.argmin(path)
                                path = tuple(np.r_[path[i_min:],path[:i_min]])
                                # to aline the orientation
                                if path[-1]<path[1]:
                                    path = path[::-1]
                                    path = tuple(np.r_[path[-1],path[:-1]])
                                if path in set_rings:
                                    pass
                                else:
                                    set_rings.add(path)
                                    
        #progress(int((ic+1)/num*100))
    #finish()
    
    return set_rings

def enumerate_king_ring(atoms_extracted, atoms_all, chemical_bond_index, flag_primitive):
    """enumerate King's rings and their primitive rings
    """

    # enumerate King's ring
    # a ring as the shortest path between two of the nearest neighbors of a given node
    set_rings = set()
    if has_igraph:        
        G = ig.Graph()
        G.add_vertices(len(atoms_all))
        G.add_edges(chemical_bond_index)
    else:
        G = nx.Graph()
        G.add_nodes_from(atoms_all)
        G.add_edges_from(chemical_bond_index)
    #progress(0)
    num = len(atoms_extracted)
    for ic, n in enumerate(atoms_extracted):
        neighbors = list(nx.all_neighbors(G, n))
        if len(neighbors)>=2:
            for n0 in neighbors:
                if has_igraph:
                    ei = G.get_eid(n, n0)
                    G.delete_edges(ei)
                else:
                    G.remove_edge(n, n0)
            combinations = list(itertools.combinations(neighbors,2))
            for comb in combinations:
                n0 = comb[0]
                n1 = comb[1]
                try:
                    if has_igraph:
                        paths = G.get_all_shortest_paths(n0, n1)
                    else:
                        paths = nx.all_shortest_paths(G, source=n0, target=n1)
                    for path in paths:
                        path.append(n)
                        path = np.array(path)
                        i_min = np.argmin(path)
                        path = tuple(np.r_[path[i_min:],path[:i_min]])
                        # to aline the orientation
                        if path[-1]<path[1]:
                            path = path[::-1]
                            path = tuple(np.r_[path[-1],path[:-1]])
                        set_rings.add(path)
                except nx.NetworkXNoPath:
                    pass #Nothing to do
            for n0 in neighbors:
                if has_igraph:
                    G.add_edges([(n, n0)])
                else:
                    G.add_edge(n, n0)
        #progress(int((ic+1)/num*100))
    
    if not(flag_primitive):
        #finish()
        return set_rings
    else:
        # remove NON-primitive rings, which have at least one shortcut
        set_primitive_rings = set()
        for path in set_rings:
            G_ring=nx.Graph()
            nx.add_cycle(G_ring, path)
            combinations = list(itertools.combinations(path,2))
            flag_primitive = True
            for comb in combinations:
                n0, n1 = comb[0], comb[1]
                l_ring = nx.shortest_path_length(G_ring, source=n0, target=n1)
                l_ori = nx.shortest_path_length(G, source=n0, target=n1)
                if l_ori<l_ring:
                    flag_primitive = False
                    break
            if flag_primitive:
                set_primitive_rings.add(tuple(path))
        #finish()
        return set_primitive_rings

def enumerate_guttman_ring(atoms_extracted, atoms_all, chemical_bond_index):
    """enumerate Guttman's rings and their primitive rings
    """    
    
    # enumerate Guttman's ring
    # a ring as the shortest path between two nearest neighbor nodes
    set_rings = set()
    if has_igraph:
        G = ig.Graph()
        G.add_vertices(len(atoms_all))
        G.add_edges(chemical_bond_index)
    else:
        G = nx.Graph()
        G.add_nodes_from(atoms_all)
        G.add_edges_from(chemical_bond_index)    
    
    bond_index = chemical_bond_index.tolist()
    
    #progress(0)
    num = chemical_bond_index.shape[0]
    for i in range(chemical_bond_index.shape[0]):
        n0 = chemical_bond_index[i,0]
        n1 = chemical_bond_index[i,1]
        if (n0 in atoms_extracted)|(n1 in atoms_extracted):            
            if has_igraph:
                ei = G.get_eid(n0, n1)
                G.delete_edges(ei)
            else:
                G.remove_edge(n0, n1)
            try:
                if has_igraph:
                    paths = G.get_all_shortest_paths(n0, n1)
                else:
                    paths = nx.all_shortest_paths(G, source=n0, target=n1)                    
                for path in paths:
                    path = np.array(path)
                    i_min = np.argmin(path)
                    path = tuple(np.r_[path[i_min:],path[:i_min]])
                    # to aline the orientation
                    if path[-1] < path[1]:
                        path = path[::-1]
                        path = tuple(np.r_[path[-1],path[:-1]])

                    set_rings.add(path)
                                    
            except nx.NetworkXNoPath:
                pass #Nothing to do
            if has_igraph:
                G.add_edges([(n0, n1)])
            else:
                G.add_edge(n0, n1)            
        #progress(int((i+1)/num*100))
        
    return set_rings

def random_rotation_matrix():
    """
        Generate ramdom rotation matrix 
    """

    x0 = np.random.random()
    y1 = 2*np.pi*np.random.random()
    y2 = 2*np.pi*np.random.random()
    r1 = np.sqrt(1.0-x0)
    r2 = np.sqrt(x0)
    u0 = np.cos(y2)*r2
    u1 = np.sin(y1)*r1
    u2 = np.cos(y1)*r1
    u3 = np.sin(y2)*r2
    coefi = 2.0*u0*u0-1.0
    coefuu = 2.0
    coefe = 2.0*u0

    R = np.zeros(shape=(3, 3))
    R[0, 0] = coefi+coefuu*u1*u1
    R[1, 1] = coefi+coefuu*u2*u2
    R[2, 2] = coefi+coefuu*u3*u3
    R[1, 2] = coefuu*u2*u3-coefe*u1
    R[2, 0] = coefuu*u3*u1-coefe*u2
    R[0, 1] = coefuu*u1*u2-coefe*u3
    R[2, 1] = coefuu*u3*u2+coefe*u1
    R[0, 2] = coefuu*u1*u3+coefe*u2
    R[1, 0] = coefuu*u2*u1+coefe*u3
    return R

class Molecule:
    """Load information of an atomic configuration
    
    The information of a configuration is loaded from a xyz or cfg (lammps) file
    
    Parameters
    ----------
    filename : str
        File name of structure model (atomic symbols and coordinates).
        
    lattice_size : float, default=False, optional 
        The length of a simulation box. Cubic is assumed in this class.
        If lattice_size is None, a non-periodic structure is assumed.

    
    Attributes
    ----------
    num_atoms : int
        The number of atoms in the lattice
        
    set_atoms : array (str)
        The set of atomic symbols included in the structure
        
    atom_symbols : array (str) of shape (num_atoms,)
        List of atomic symbols
    
    xyz : array (float) of shape (num_atoms, 3)
        Atomic (x,y,z) coordinates
        
    chemical_bond_index_atoms : array (int) of shape (# of bonds, 2)
        List of atom pairs which have a chemical bonds
    
    chemical_bond_index_cells : array (int) of shape(# of bonds, 3)
        Indices of a super cell that the second atom of a bond is included.
        The element must be one of {-1, 0, +1}.
    
    bond_pair_atom_symbols : list of list (str) of shape (2,)
        List of paird atom symbols to make chemical bonds. 
        
    bond_pair_dist_max : list of list (int)
        List of maximum lengths of chemical bonds between an atom pair. 
        
    bond_flag_periodicity : Boolean 
        Whether periodic condition is assumed.    
    """

    def __init__(self, filename=None, lattice_size=False):
        
        self.unit_atoms_index = None
        self.ghost = None
        
        if filename is None:
            self.results = None
        else:            
            _, ext = os.path.splitext(filename)
            if ext=='.xyz':
                self.read_xyz(filename)
                self.lattice_size = lattice_size
            elif ext=='.cfg':
                self.read_cfg_lammps(filename)
                self.lattice_size = lattice_size
            elif ext=='.npz':
                data = np.load(filename, allow_pickle=True)
                self.atom_symbols = data['atom_symbols']
                self.xyz = data['xyz']
                self.num_atoms = self.xyz.shape[0]
                self.lattice_size = data['lattice_size']
                
                if ('chemical_bond_index_atoms' in data.files) and (data['chemical_bond_index_atoms'] is not None):
                    self.chemical_bond_index_atoms = data['chemical_bond_index_atoms']
                if ('chemical_bond_index_cells' in data.files) and (data['chemical_bond_index_cells'] is not None):
                    self.chemical_bond_index_cells = data['chemical_bond_index_cells']
                if ('bond_pair_atom_symbols' in data.files) and (data['bond_pair_atom_symbols'] is not None):
                    self.bond_pair_atom_symbols = data['bond_pair_atom_symbols']
                if ('bond_pair_dist_max' in data.files) and (data['bond_pair_dist_max'] is not None):
                    self.bond_pair_dist_max = data['bond_pair_dist_max']
                if ('bond_flag_periodicity' in data.files)  and (data['bond_flag_periodicity'] is not None):
                    self.bond_flag_periodicity = data['bond_flag_periodicity']            
    
    def set_results(self, results, lattice_size=False):
        """
        set results data.

        Parameters
        ----------
        results : class
            core.data.Results class object.
        lattice_size : float, optional
            lattice size. The default is False.

        Returns
        -------
        None.

        """
        self.results = results
        
        if results is None:
            return
        
        self.num_atoms = self.results.atoms.number
        self.xyz = self.results.atoms.positions
        self.atom_symbols = np.array(self.results.atoms.elements)
        self.set_atoms = np.unique(self.atom_symbols)
        if lattice_size is False:
            # TODO temporary for Cubic system
            self.lattice_size = self.results.atoms.volume.vectors[0][0]*2.
        else:         
            self.lattice_size = lattice_size
        self.xyz = self.xyz % self.lattice_size
    
    def set_atoms(self, atoms, lattice_size=False):        
        self.atoms = atoms        
        if atoms is None:
            return
        
        self.num_atoms = self.atoms.number
        self.xyz = self.atoms.positions
        self.atom_symbols = np.array(self.atoms.elements)
        #self.set_atoms = np.unique(self.atom_symbols)
        """
        if lattice_size is False:
            # TODO temporary for Cubic system
            self.lattice_size = self.atoms.volume.vectors[0][0]*2.
        else:         
            self.lattice_size = lattice_size
        self.xyz = self.xyz % self.lattice_size
        """
        
    def read_xyz(self, filename):
        """Load structure information from a xyz file
        
        Load structure information of atomic symbols and coordinates from a xyz file.
        
        Parameters
        ----------
        filename : str
            xyz file name.
        """

        with open(filename,'r') as f:
            self.num_atoms = int(f.readline().rstrip()) #the number of atoms in the structure
            f.readline()
            atom_symbols = list() #list of atomic symbols
            self.xyz = np.zeros((self.num_atoms,3))
            for n in range(self.num_atoms):
                line = f.readline().rstrip()
                line = line.split()
                atom_symbols.append(line[0]) #atomic symbol
                self.xyz[n,:] = np.array(line[1:4]).astype(np.float64) #atomic coordinate (x,y,z)
            self.atom_symbols = np.array(atom_symbols)
            self.set_atoms = np.unique(self.atom_symbols) #the set of atomic symbols

    def save_to_npz(self, filename):
        """Save structure information

        Save the atomic configuration and chemical bonds.

        Parameters
        ----------
        filename : str
            Saved file name.
        """
        
        if not(hasattr(self,'chemical_bond_index_atoms')):
            self.chemical_bond_index_atoms = None
        if not(hasattr(self,'chemical_bond_index_cells')):
            self.chemical_bond_index_cells = None
        if not(hasattr(self,'bond_pair_atom_symbols')):
            self.bond_pair_atom_symbols = None
        if not(hasattr(self,'bond_pair_dist_max')):
            self.bond_pair_dist_max = None
        if not(hasattr(self,'bond_flag_periodicity')):
            self.bond_flag_periodicity = None
        
        np.savez_compressed(filename,\
            atom_symbols = self.atom_symbols,\
            xyz = self.xyz,\
            lattice_size = self.lattice_size, \
            chemical_bond_index_atoms = self.chemical_bond_index_atoms,\
            chemical_bond_index_cells = self.chemical_bond_index_cells,\
            bond_pair_atom_symbols = self.bond_pair_atom_symbols,\
            bond_pair_dist_max = self.bond_pair_dist_max,\
            bond_flag_periodicity = self.bond_flag_periodicity
            )

    def save_to_xyz(self, filename, index_atoms=None):
        """Save the atomic configuration to a xyz file
        
        Save the configuration of all or a subset of atoms to a xyz file.

        Parameters
        ----------
        filename : str
            xyz file name to save.
            
        index_atoms : array (int), default=None
            Indices of atoms to save.
        """
        
        #If index_atoms is None, save all atoms
        if index_atoms is None:
            index_atoms = np.arange(self.num_atoms)
        else:
            index_atoms = np.array(index_atoms)
        N = index_atoms.size
        
        #Save the atomic configuration to a xyz file 
        with open(filename, 'w') as f:
            f.write(str(N)+'\n\n')
            N = len(index_atoms)
            for n in range(N):
                id = index_atoms[n]
                x = self.xyz[id,0]
                y = self.xyz[id,1]
                z = self.xyz[id,2]
                line = '{:} {:.3f} {:.3f} {:.3f}'.format(self.atom_symbols[id], x, y, z)
                f.write(line+'\n')

    def make_bond_periodic(self, pair_atom_symbols, pair_dist_max, flag_periodicity=False, p_pair=0.3):
        
        dist_max = np.max(pair_dist_max)
        grid = self.atoms.grid
        
        #progress(0)

        index_atoms = []
        index_cells = []
        distance = []
        for ic in range(self.atoms.number):        
            grid.neighbours(ic, dist_max, 0, self.atoms.number-1)
            for inei, d, shift in zip(grid.inei, grid.d, np.array(grid.shifts).T):
                if ic < inei:
                    index_atoms.append([ic, inei])
                    index_cells.append([shift[0],shift[1],shift[2]])
                    distance.append(d)
            
            #progress(int((ic+1)/self.results.atoms.number*100))
        
        index_atoms = np.array(index_atoms)
        index_cells = np.array(index_cells)
        distance = np.array(distance)
        
        if len(index_atoms) == 0:
            print('Warning, not found index atoms : molecule.make_bond_periodic')
            return
                
        # judge if the distance of each pairs of atoms is less than pair_dist
        bool_bond = np.zeros(distance.shape, dtype=bool)
        for (a0, a1), dist_max in zip(pair_atom_symbols, pair_dist_max):
            bool_bond = (((self.atom_symbols[index_atoms[:,0]]==a0) & (self.atom_symbols[index_atoms[:,1]]==a1)) \
                | ((self.atom_symbols[index_atoms[:,0]]==a1) & (self.atom_symbols[index_atoms[:,1]]==a0))) \
                & (distance < dist_max) | bool_bond
        self.chemical_bond_index_atoms = index_atoms[bool_bond, :]
        if flag_periodicity:
            self.chemical_bond_index_cells = index_cells[bool_bond, :]
        
        #save options to make bonds
        self.bond_pair_atom_symbols = pair_atom_symbols 
        self.bond_pair_dist_max = pair_dist_max
        self.bond_flag_periodicity = flag_periodicity
        
    
    def make_bond_ghost(self,pair_atom_symbols, pair_dist_max, flag_periodicity=False, 
                        p_pair=0.3, ghost_cell=1):

        dist_max = np.max(pair_dist_max)
        
        #progress(0)
        
        # make supercell position 3x3x3
        _positions = []
        self.ghost = []
        for i, x in enumerate(self.atoms.norm_positions):
            _positions.append(x)
            self.ghost.append(i)
        
        for ix in [-ghost_cell,0,ghost_cell]:
            for iy in [-ghost_cell,0,ghost_cell]:
                for iz in [-ghost_cell,0,ghost_cell]:
                    if ix == 0 and iy == 0 and iz == 0:
                        continue
                    for i, x in enumerate(self.atoms.norm_positions):
                        pos = x + np.array([ix*2,iy*2,iz*2])
                        _positions.append(pos)
                        self.ghost.append(i)
        
        nghost = 2*ghost_cell+1
        positions = np.array(_positions)
        positions = positions + np.array([float(nghost),float(nghost),float(nghost)])
        positions /= float(nghost)
        positions -= 1.
        vectors = np.array(self.atoms.volume.vectors)
        self.num_atoms = positions.shape[0]
        
        m = np.identity(3)*nghost
        vectors *= m
        grid = Grid(positions, vectors)

        self.unit_atoms_index = [i for i in range(self.atoms.number)]

        index_atoms = []
        index_cells = []
        distance = []
        n = len(positions)
        for ic in range(n):        
            grid.neighbours(ic, dist_max, 0, n-1)
            for inei, d, shift in zip(grid.inei, grid.d, np.array(grid.shifts).T):
                if ic < inei:
                    index_atoms.append([ic, inei])
                    index_cells.append([shift[0],shift[1],shift[2]])
                    distance.append(d)
        
        type_atoms = [[self.ghost[i], self.ghost[j]] for i, j in index_atoms]
        
        #print(index_atoms)
            #progress(int((ic+1)/self.results.atoms.number*100))

        type_atoms = np.array(type_atoms)
        index_atoms = np.array(index_atoms)
        index_cells = np.array(index_cells)
        distance = np.array(distance)
                
        # judge if the distance of each pairs of atoms is less than pair_dist
        bool_bond = np.zeros(distance.shape, dtype=bool)
        for (a0, a1), dist_max in zip(pair_atom_symbols, pair_dist_max):
            bool_bond = (((self.atom_symbols[type_atoms[:,0]]==a0) & (self.atom_symbols[type_atoms[:,1]]==a1)) \
                | ((self.atom_symbols[type_atoms[:,0]]==a1) & (self.atom_symbols[type_atoms[:,1]]==a0))) \
                & (distance < dist_max) | bool_bond
        self.chemical_bond_index_atoms = index_atoms[bool_bond, :]
        if flag_periodicity:
            self.chemical_bond_index_cells = index_cells[bool_bond, :]
        
        #save options to make bonds
        self.bond_pair_atom_symbols = pair_atom_symbols 
        self.bond_pair_dist_max = pair_dist_max
        self.bond_flag_periodicity = flag_periodicity
    
    def make_bond(self, pair_atom_symbols, pair_dist_max, flag_periodicity=False, p_pair=0.3):
        """Enumerate chemical bonds
            
            Enumerate all chemical bonds between atom pairs indicated by pair_atom_symbols
            and whose distance is less than pair_dist_max. If flag_periodicity is True,
            the periodic condition is assumed and then chemical bonds over the neigrbor 
            cells (lattiecs) are generated.

        Parameters
        ----------
        pair_atom_symbols : list of list (str) of shape (2,)
            List of paird atom symbols to make chemical bonds.
        pair_dist_max : list of list (int)
            List of maximum lengths of chemical bonds between an atom pair.
        flag_periodicity : Boolean, default False
            Whether periodic condition is assumed.
        p_pair : float
            The ratio of the number of atom pairs.
            p_pair is used to make a variable with fixed size in jit_distance.
        """

        # enumerate all pairs of atoms whose distance is less than dist_max 
        dist_max = np.max(pair_dist_max)
        if flag_periodicity:
            if self.lattice_size is False:
                print('Molecule.lattice_size is not indicated')
            else:
                # distance, index_atoms, index_cells = jit_distance_supercell(self.xyz.astype(np.float32), self.lattice_size, dist_max.astype(np.float32), np.float32(p_pair))
                distance, index_atoms, index_cells = distance_supercell(self.xyz, self.lattice_size, dist_max, p_pair)
        else:
            distance, index_atoms = jit_distance(self.xyz.astype(np.float32), dist_max.astype(np.float32), np.float32(p_pair))
                
        # judge if the distance of each pairs of atoms is less than pair_dist
        bool_bond = np.zeros(distance.shape, dtype=bool)
        for (a0, a1), dist_max in zip(pair_atom_symbols, pair_dist_max):
            bool_bond = (((self.atom_symbols[index_atoms[:,0]]==a0) & (self.atom_symbols[index_atoms[:,1]]==a1)) \
                | ((self.atom_symbols[index_atoms[:,0]]==a1) & (self.atom_symbols[index_atoms[:,1]]==a0))) \
                & (distance < dist_max) | bool_bond
        self.chemical_bond_index_atoms = index_atoms[bool_bond, :]
        if flag_periodicity:
            self.chemical_bond_index_cells = index_cells[bool_bond, :]
            
        #save options to make bonds
        self.bond_pair_atom_symbols = pair_atom_symbols 
        self.bond_pair_dist_max = pair_dist_max
        self.bond_flag_periodicity = flag_periodicity
        
class RING(Molecule):
    
    class RingType:
        GUTTMAN        = 0 # Guttman
        KING           = 1 # King
        PRIMITIVE      = 2 # Primitive
        PRIMITIVE_KING = 3 # Primitive-King
    
    """Analysis of the topology of chemical bonds
    
    This class can be used for topological analysis of network forming chemical structure
    based on rings and tetrahedra.
    
    Parameters
    ----------
    filename : str
        File name of structure model (atomic symbols and coordinates).
        
    lattice_size : float, default=None, optional 
        The length of a simulation box. Cubic is assumed in this class.
        If lattice_size is None, a non-periodic structure is assumed.

    
    Attributes
    ----------
    num_atoms : int
        The number of atoms in the lattice
        
    set_atoms : array (str)
        The set of atomic symbols included in the structure
        
    atom_symbols : array (str) of shape (num_atoms,)
        List of atomic symbols
    
    xyz : array (float) of shape (num_atoms, 3)
        Atomic (x,y,z) coordinates
        
    chemical_bond_index_atoms : array (int) of shape (# of bonds, 2)
        List of atom pairs which have a chemical bonds
    
    chemical_bond_index_cells : array (int) of shape(# of bonds, 3)
        Indices of a super cell that the second atom of a bond is included.
        The element must be one of {-1, 0, +1}.
    
    bond_pair_atom_symbols : list of list (str) of shape (2,)
        List of paird atom symbols to make chemical bonds. 
        
    bond_pair_dist_max : list of list (int)
        List of maximum lengths of chemical bonds between an atom pair. 
        
    bond_flag_periodicity : Boolean 
        Whether periodic condition is assumed.
        
    ring_type : str, {'Primitive', 'King', 'Primitive_King', 'Guttman'}
        A type of enumerated rings.
    
    rings : list of list (int) of shape (# of atoms in a ring)
        List of the sets of atom index in rings
    
    ring_centers : array of shape (# of rings, 3)
        Ring center coordinates computate from the average over atoms in a ring.
    
    ring_normal_vectors : array of shape (# of rings, 3)
        Normal vectors of planes fitted by atoms in each ring.
        
    ring_ellipse_long_vectors :  array of shape (# of rings, 3)
        Long axes of ellipses fitted by atomic configuration in each ring.
    
    ring_ellipse_short_vectors :  array of shape (# of rings, 3)
        Long axes of ellipses fitted by atomic configuration in each ring.
        
    ring_ellipsoid_lengths : array of shape (# of rings, 3)
        Length of three axes of an ellipsoid fitted by atomic configuration in a ring.
        
    ring_near_pairs : array of shape (# of ring pairs, 2)
        Pairs of ring indeces whose distance is less than a threshold.
    
    ring_pair_distances : array of shape (# of ring pairs,)
        Distances of ring centers.
    
    ring_pair_index_cells : array of shape (# of ring pairs, 3)
        Indices of a super cell that the second ring center is included.
        The element must be one of {-1, 0, +1}.
    
    ring_circleness : array of shape (# of rings,)
        Circleness of rings.
    
    tetra_q_values : array (int) of shape (# of tetrahedra,)
        q-values of tetrahedrons, which evaluate symmetry.
    
    tetra_center_indices : array (int) of shape (# of tetrahedra,)
        Indices of center atoms of tetrahedra.
    
    tetra_neighbor_indices : array (int) of shape (# of tetrahedra, 4)
        Indices of vertex atoms of tetrahedra.
    
    tetra_neighbor_distances : array (int) of shape (# of tetrahedra, 4)
        Distances between a center atom and a tetrahedral vertex.

    tetra_neighbor_cell_indices : array (int) of shape (# of tetrahedra, 4, 3)
        Cell indices {-1, 0, +1}, from center atoms, of a tetrahedral vertex.
    """
    
    def __init__(self, filename=None, lattice_size=False):
        
        if filename is None:
            super().__init__()
        else:
            _, ext = os.path.splitext(filename)
            if ext=='.npz':
                data = np.load(filename, allow_pickle=True)
                
                #basic data (lattice and atomic configuration)
                self.atom_symbols = data['atom_symbols']
                self.xyz = data['xyz']
                self.num_atoms = self.xyz.shape[0]
                self.lattice_size = data['lattice_size']
                self.set_atoms = np.array(list(set(self.atom_symbols)))
                
                #chemical bonds
                if ('chemical_bond_index_atoms' in data.files) and (data['chemical_bond_index_atoms'] is not None):
                    self.chemical_bond_index_atoms = data['chemical_bond_index_atoms']
                if ('chemical_bond_index_cells' in data.files) and (data['chemical_bond_index_cells'] is not None):
                    self.chemical_bond_index_cells = data['chemical_bond_index_cells']
                if ('bond_pair_atom_symbols' in data.files) and (data['bond_pair_atom_symbols'] is not None):
                    self.bond_pair_atom_symbols = data['bond_pair_atom_symbols']
                if ('bond_pair_dist_max' in data.files) and (data['bond_pair_dist_max'] is not None):
                    self.bond_pair_dist_max = data['bond_pair_dist_max']
                if ('bond_flag_periodicity' in data.files) and (data['bond_flag_periodicity'] is not None):
                    self.bond_flag_periodicity = data['bond_flag_periodicity']
                    
                #rings
                if ('ring_type' in data.files) and (data['ring_type'] is not None):
                    self.ring_type = data['ring_type']
                if ('rings' in data.files) and (data['rings'] is not None):
                    self.rings = data['rings']
                    
                #ring characters
                if ('ring_centers' in data.files) and (data['ring_centers'] is not None):
                    self.ring_centers = data['ring_centers']
                if ('ring_normal_vectors' in data.files) and (data['ring_normal_vectors'] is not None):
                    self.ring_normal_vectors = data['ring_normal_vectors']
                if ('ring_near_pairs' in data.files) and (data['ring_near_pairs'] is not None):
                    self.ring_near_pairs = data['ring_near_pairs']
                if ('ring_pair_distances' in data.files) and (data['ring_pair_distances'] is not None):
                    self.ring_pair_distances = data['ring_pair_distances']
                if ('ring_pair_index_cells' in data.files) and (data['ring_pair_index_cells'] is not None):
                        self.ring_pair_index_cells = data['ring_pair_index_cells']
                if ('ring_circleness' in data.files) and (data['ring_circleness'] is not None):
                        self.ring_circleness = data['ring_circleness']
                if ('ring_ellipse_long_vectors' in data.files) and (data['ring_ellipse_long_vectors'] is not None):
                        self.ring_ellipse_long_vectors = data['ring_ellipse_long_vectors']
                if ('ring_ellipse_short_vectors' in data.files) and (data['ring_ellipse_short_vectors'] is not None):
                        self.ring_ellipse_short_vectors = data['ring_ellipse_short_vectors']
                if ('ring_ellipsoid_lengths' in data.files) and (data['ring_ellipsoid_lengths'] is not None):
                        self.ring_ellipsoid_lengths = data['ring_ellipsoid_lengths']
    
                #tetrahedra
                if ('tetra_q_values' in data.files) and (data['tetra_q_values'] is not None):
                    self.tetra_q_values = data['tetra_q_values']
                if ('tetra_center_indices' in data.files) and (data['tetra_center_indices'] is not None):
                    self.tetra_center_indices = data['tetra_center_indices']
                if ('tetra_neighbor_indices' in data.files) and (data['tetra_neighbor_indices'] is not None):
                    self.tetra_neighbor_indices = data['tetra_neighbor_indices']
                if ('tetra_neighbor_distances' in data.files) and (data['tetra_neighbor_distances'] is not None):
                    self.tetra_neighbor_distances = data['tetra_neighbor_distances']
                if ('tetra_neighbor_cell_indices' in data.files) and (data['tetra_neighbor_cell_indices'] is not None):
                    self.tetra_neighbor_cell_indices = data['tetra_neighbor_cell_indices']
            elif (ext=='.xyz')|(ext=='.cfg'):
                #load data from a xyz or cfg file
                super().__init__(filename, lattice_size)
            else:
                print('An input file must be *.xyz, *.cfg or *.npz!')

    def save_to_npz(self, filename):
        """Save data of computed rings and tetrahedra
        
        Save data of computed rings and tetrahedra and their characters.
        
        Parameters
        ----------
        filename : (str)
            Name of a saved npz file. 
        """
        _, ext = os.path.splitext(filename)
        if ext=='.npz':
            #save computed data to file '.npz' file
            if not(hasattr(self,'chemical_bond_index_atoms')):
                self.chemical_bond_index_atoms = None
            if not(hasattr(self,'chemical_bond_index_cells')):
                self.chemical_bond_index_cells = None
            
            if not(hasattr(self,'bond_pair_atom_symbols')):
                self.bond_pair_atom_symbols = None
            if not(hasattr(self,'bond_pair_dist_max')):
                self.bond_pair_dist_max = None
            if not(hasattr(self,'bond_flag_periodicity')):
                self.bond_flag_periodicity = None    
            
            if not(hasattr(self, 'ring_type')):
                self.ring_type = None
            if not(hasattr(self, 'rings')):
                self.rings = None
            if not(hasattr(self, 'ring_centers')):
                self.ring_centers = None
            if not(hasattr(self, 'ring_normal_vectors')):
                self.ring_normal_vectors = None
            if not(hasattr(self, 'ring_near_pairs')):
                self.ring_near_pairs = None
            if not(hasattr(self, 'ring_pair_distances')):
                self.ring_pair_distances = None
            if not(hasattr(self, 'ring_pair_index_cells')):
                self.ring_pair_index_cells = None
            if not(hasattr(self, 'ring_circleness')):
                self.ring_circleness = None
            if not(hasattr(self, 'ring_ellipse_long_vectors')):
                self.ring_ellipse_long_vectors = None
            if not(hasattr(self, 'ring_ellipse_short_vectors')):
                self.ring_ellipse_short_vectors = None
            if not(hasattr(self, 'ring_ellipsoid_lengths')):
                self.ring_ellipsoid_lengths = None

            if not(hasattr(self, 'tetra_q_values')):
                self.tetra_q_values = None
            if not(hasattr(self, 'tetra_center_indices')):
                self.tetra_center_indices = None
            if not(hasattr(self, 'tetra_neighbor_indices')):
                self.tetra_neighbor_indices = None
            if not(hasattr(self, 'tetra_neighbor_distances')):
                self.tetra_neighbor_distances = None
            if not(hasattr(self, 'tetra_neighbor_cell_indices')):
                self.tetra_neighbor_cell_indices = None
            
            np.savez_compressed(filename,\
                atom_symbols = self.atom_symbols,\
                xyz = self.xyz,\
                lattice_size = self.lattice_size, \
                chemical_bond_index_atoms = self.chemical_bond_index_atoms,\
                chemical_bond_index_cells = self.chemical_bond_index_cells,\
                bond_pair_atom_symbols = self.bond_pair_atom_symbols,\
                bond_pair_dist_max = self.bond_pair_dist_max,\
                bond_flag_periodicity = self.bond_flag_periodicity,\
                ring_type = self.ring_type,\
                rings = self.rings,\
                ring_centers = self.ring_centers,\
                ring_normal_vectors = self.ring_normal_vectors,\
                ring_near_pairs = self.ring_near_pairs,\
                ring_pair_distances = self.ring_pair_distances, \
                ring_pair_index_cells = self.ring_pair_index_cells, \
                ring_circleness = self.ring_circleness, \
                ring_ellipse_long_vectors = self.ring_ellipse_long_vectors, \
                ring_ellipse_short_vectors = self.ring_ellipse_short_vectors, \
                ring_ellipsoid_lengths = self.ring_ellipsoid_lengths, \
                tetra_q_values = self.tetra_q_values, \
                tetra_center_indices = self.tetra_center_indices, \
                tetra_neighbor_indices = self.tetra_neighbor_indices, \
                tetra_neighbor_distances = self.tetra_neighbor_distances, \
                tetra_neighbor_cell_indices = self.tetra_neighbor_cell_indices)
        else:
            print("Error: The file extension must be \".npz\"!")
            
    def calculate(self, ring_type, pair_atom_symbols, pair_dist_max, periodicity=True, p_pair=0.3,
                  cutoff_size=None, atoms_extracted=None, chain=False):
        
        num_parallel = 0
        if periodicity:            
            self.make_bond_periodic(pair_atom_symbols, pair_dist_max, 
                                    periodicity, p_pair)
            self.enumerate_ring(ring_type, cutoff_size, num_parallel, atoms_extracted)
                        
            if chain:
                rings = self.rings
            else:                
                rings = self.select_periodic_rings()
                if len(rings) == 0:
                    self.make_bond_ghost(pair_atom_symbols, pair_dist_max, 
                                            periodicity, p_pair)
                    self.enumerate_ring(ring_type, cutoff_size, num_parallel, atoms_extracted)
                    rings = self.select_periodic_rings()
        else:
            self.make_bond(pair_atom_symbols, pair_dist_max)
            self.enumerate_ring(ring_type, cutoff_size, num_parallel, atoms_extracted)
            rings = self.rings
            
        self.rings = [list(r) for r in list(rings)] # save rings as list
        
        return self.rings
            
    def enumerate_ring(self, ring_type, cutoff_size=None, num_parallel=0, atoms_extracted=None):
        """Enumerate rings
        
        Enumerate rings from a network generated from chemical bonds
        
        Parameters
        ----------
        ring_type : str, {'Primitive', 'King', 'Primitive_King', 'Guttman'}
            Type of enumerated rings.
            Choose one from Primitive, King, Primitive_King, and Guttman.
            
        cutoff_size : int, optional
            The maximum size of enumerated rings, by default None
            This is a parameter only for primitive ring enumeration. 
        
        num_parallel : int, optional
            The number of CPU cores to implement the enumeration, by default 0
            For non parallel computation, set 0.
            To use the maximum numbers of CPU cores, set -1. 
            The number of CPU cores in your computer can be checked by 'os.cpu_count()'.
            
        atoms_extracted : array (int), optional
            Subset of atoms used as start nodes of ring enumeration.
            This option should not be used (set to None) if the periodical condition is assumed.
            This option is used only when bond_flag_periodicity=True. 
        """

        atoms_all = np.arange(self.num_atoms)        
        if atoms_extracted is None:
            atoms_extracted = atoms_all
        if self.unit_atoms_index is not None:            
            atoms_extracted = self.unit_atoms_index
                
        flag_type = True
        if (ring_type == RING.RingType.PRIMITIVE)&(num_parallel==0):
            set_rings = enumerate_primitive_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, cutoff_size)
        # elif (ring_type=='Primitive')&(num_parallel!=0):
        #     set_rings = parallel_enumerate_primitive_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms,\
        #         cutoff_size, num_parallel=num_parallel)
        elif (ring_type == RING.RingType.PRIMITIVE_KING)&(num_parallel==0):
            set_rings = enumerate_king_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, flag_primitive=True)
        # elif (ring_type=='Primitive_King')&(num_parallel!=0):
        #     set_rings = parallel_enumerate_king_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms,\
        #         flag_primitive=True, num_parallel=num_parallel)
        elif (ring_type == RING.RingType.KING)&(num_parallel==0):
            set_rings = enumerate_king_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, flag_primitive=False)
        # elif (ring_type=='King')&(num_parallel!=0):
        #     set_rings = parallel_enumerate_king_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, \
        #         flag_primitive=False, num_parallel=num_parallel)
        elif (ring_type == RING.RingType.GUTTMAN)&(num_parallel==0):
            set_rings = enumerate_guttman_ring(atoms_extracted, atoms_all,self.chemical_bond_index_atoms)
        # elif (ring_type=='Guttman')&(num_parallel!=0):
        #     set_rings = parallel_enumerate_guttman_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, \
        #         num_parallel=num_parallel)
        else:
            print('The indicated ring type is incorrect.')
            print('Choose a fring tyep from Primitive, King, Primitive_King, and Guttman.') 
            flag_type=False

        if flag_type:
            self.ring_type = ring_type
            self.rings = [list(r) for r in list(set_rings)] # save rings as list
    
    def select_periodic_rings(self):
        _rings = set()
        
        # check summation shift vector
        bond_index = self.chemical_bond_index_atoms.tolist()
        for i, ring in enumerate(self.rings):
            sum_shift = np.zeros(3, dtype='int8')
            lpath = list(ring)
            for ns, ne in zip(lpath, lpath[1:]+lpath[:1]):
                sign = 1
                if ns > ne:
                    ns, ne = ne, ns
                    sign = -1    
                ib = bond_index.index([ns,ne])
                shift = self.chemical_bond_index_cells[ib]
                sum_shift += sign*shift            
            if not np.all(sum_shift == 0):
                continue
            
            if self.ghost is not None:
                path = [self.ghost[i] for i in ring] # create ghost path
            else:
                path = ring                
            path = np.array(path)
            i_min = np.argmin(path)
            path = tuple(np.r_[path[i_min:],path[:i_min]])
            # to aline the orientation
            if path[-1] < path[1]:
                path = path[::-1]
                path = tuple(np.r_[path[-1],path[:-1]])
            _rings.add(path)
            
        return _rings
    
    """
    def check_bond(self):        
        self.shifts = []
        for ring in self.rings:            
            ring_list = np.array(ring)
            shifts = []        
            centres = self.atoms.norm_positions        
            for i, j in zip(ring_list.tolist(), np.roll(ring_list,-1).tolist()):
                shift = centres[j]-centres[i]+3.            
                shift = [int(x/2.)-1 for x in shift]            
                shifts.append(shift)
            self.shifts.append(shifts)
    """        