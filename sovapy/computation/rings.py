import numpy as np
from numpy import linalg
import time, itertools, os, shutil, pickle, multiprocessing
from concurrent.futures import ThreadPoolExecutor
from concurrent import futures
import networkx as nx
from ..core.gridding import Grid
from functools import partial
from tqdm import tqdm
import copy
#from tqdm.notebook import tqdm

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

def sub_parallel_enumerate_primitive_ring(n_source, G, D, cutoff_size):
    set_rings=set()
        
    # print(n_source,flush=True)
    length = D[n_source]
    for L in range(1,round(cutoff_size/2)+1):
        i = np.array(list(length.values()))==L
        node_check = np.array(list(length.keys()))[i]
        
        # rings with even number of nodes
        for n_target in node_check:
            paths = [p for p in nx.all_shortest_paths(G,source=n_source,target=n_target)]
            if (len(paths)>1):
                #enumerate pairs of paths
                for p0, p1 in itertools.combinations(paths,2):
                    if len(set(p0)&set(p1))==2: #common nodes are only n_source and n_target

                        flag = True
                        set_node = p0 + p1[-2:0:-1]
                        num_ring_tmp = len(set_node)
                        dic_id_v = dict()
                        for n_tmp, v_tmp in enumerate(set_node):
                            dic_id_v[v_tmp] = n_tmp

                        for ni,nj in  itertools.combinations(set_node,2):
                            i_ni, i_nj = dic_id_v[ni], dic_id_v[nj]
                            if i_ni<i_nj:
                                dist_ring_tmp = min([i_nj - i_ni, i_ni + num_ring_tmp - i_nj]) 
                            else:
                                dist_ring_tmp = min([i_ni - i_nj, i_nj + num_ring_tmp - i_ni])
                            if D[ni][nj] < dist_ring_tmp:
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
                            set_rings.add(path)

        # rings with odd number of nodes
        for nt0, nt1 in itertools.combinations(node_check,2):
            if nt0 in G.neighbors(nt1): # if nt0 and nt1 are linked
                paths0 = [p for p in nx.all_shortest_paths(G,source=n_source,target=nt0)]
                paths1 = [p for p in nx.all_shortest_paths(G,source=n_source,target=nt1)]
                for p0,p1 in itertools.product(paths0, paths1):
                    if len(set(p0)&set(p1))==1: #common nodes are only n_source and n_target

                        flag = True
                        set_node = p0 + p1[-2:0:-1]
                        num_ring_tmp = len(set_node)
                        dic_id_v = dict()
                        for n_tmp, v_tmp in enumerate(set_node):
                            dic_id_v[v_tmp] = n_tmp

                        for ni,nj in  itertools.combinations(set_node,2):
                            i_ni, i_nj = dic_id_v[ni], dic_id_v[nj]
                            if i_ni<i_nj:
                                dist_ring_tmp = min([i_nj - i_ni, i_ni + num_ring_tmp - i_nj]) 
                            else:
                                dist_ring_tmp = min([i_ni - i_nj, i_nj + num_ring_tmp - i_ni])
                            if D[ni][nj] < dist_ring_tmp:
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
    return set_rings

def parallel_enumerate_primitive_ring(atoms_extracted, atoms_all, chemical_bond_index, cutoff_size, num_parallel=-1):
    """Enumerate primitive rings
    
    Parameters
    ----------
    num_atoms : int
        the number of atoms
    chemical_bond_index : list
        the list of atom pairs with bonded
    cutoff_size : int
        the maximum number of atoms in enumerated rings
    num_parallel : int, optional
        the number of parallel numbers, by default -1
        When -1, the maximum number of cores (os.cpu_count-2) is set. 
    
    Returns
    -------
    set_rings_all: set (list of int)
        the enumerated primitive rings
    """

    # the number of processes in parallel computation
    if num_parallel > os.cpu_count():
        num_parallel = os.cpu_count()-2
    elif num_parallel == -1:
        num_parallel = os.cpu_count()-2

    list_nodes = np.array_split(atoms_extracted, num_parallel)
    
    # initialize graph object
    G=nx.Graph()
    G.add_nodes_from(atoms_all)
    G.add_edges_from(chemical_bond_index)
    D = dict(nx.all_pairs_shortest_path_length(G, cutoff=round(cutoff_size/2)+1))

    # make function with fixed parameters
    fun_enum_rings = partial(sub_parallel_enumerate_primitive_ring, G=G, D=D, cutoff_size=cutoff_size) 
    
    num = len(atoms_extracted)
    
    with tqdm(total=num) as progress_bar:
        
        with ThreadPoolExecutor() as executor:
            
            set_rings_all = set()
            procs = []
            for n in atoms_extracted:
                proc = executor.submit(fun_enum_rings, n)            
                procs.append(proc)
                
            for f in futures.as_completed(procs):
                progress_bar.update(1)
                                
            executor.shutdown(wait=True)
            
            for p in procs:
                s= p.result()            
                if len(s)>0:
                    for ss in list(s):
                        set_rings_all.add(ss)
    
    return set_rings_all

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
    progress_bar = tqdm(total = num)
    
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
        progress_bar.update()
        
    #finish()
    
    return set_rings

def sub_parallel_enumerate_king_ring(n, G, flag_primitive):
    set_rings = set()    
    _G = copy.deepcopy(G)
    
    neighbors = list(nx.all_neighbors(_G, n))
    if len(neighbors)>=2:
        for n0 in neighbors:
            if has_igraph:
                ei = _G.get_eid(n, n0)
                _G.delete_edges(ei)
            else:
                _G.remove_edge(n, n0)
        combinations = list(itertools.combinations(neighbors,2))
        for comb in combinations:
            n0 = comb[0]
            n1 = comb[1]
            try:
                if has_igraph:
                    paths = _G.get_all_shortest_paths(n0, n1)
                else:
                    paths = nx.all_shortest_paths(_G, source=n0, target=n1)
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
                _G.add_edges([(n, n0)])
            else:
                _G.add_edge(n, n0)        
    
    if not(flag_primitive):        
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
                l_ori = nx.shortest_path_length(_G, source=n0, target=n1)
                if l_ori<l_ring:
                    flag_primitive = False
                    break
            if flag_primitive:
                set_primitive_rings.add(tuple(path))        
        return set_primitive_rings

def parallel_enumerate_king_ring(atoms_extracted, atoms_all, chemical_bond_index, flag_primitive, num_parallel=-1):

    # the number of processes in parallel computation
    if num_parallel > os.cpu_count():
        num_parallel = os.cpu_count()-2
    elif num_parallel == -1:
        num_parallel = os.cpu_count()-2

    #list_nodes = np.array_split(atoms_extracted, num_parallel)
    num = len(atoms_extracted)

    # initialize graph object
    if has_igraph:        
        G = ig.Graph()
        G.add_vertices(len(atoms_all))
        G.add_edges(chemical_bond_index)
    else:
        G = nx.Graph()
        G.add_nodes_from(atoms_all)
        G.add_edges_from(chemical_bond_index)

    # make function with fixed parameters
    fun_enum_rings = partial(sub_parallel_enumerate_king_ring, G=G, flag_primitive=flag_primitive) 
    
    with tqdm(total=num) as progress_bar:
        
        with ThreadPoolExecutor() as executor:
            #progress_bar = tqdm(total = num)
            set_rings_all = set()
            procs = []
            for n in atoms_extracted:
                proc = executor.submit(fun_enum_rings, n)            
                procs.append(proc)
                
            for f in futures.as_completed(procs):
                progress_bar.update(1)
                                
            executor.shutdown(wait=True)
            
            for p in procs:
                s= p.result()            
                if len(s)>0:
                    for ss in list(s):
                        set_rings_all.add(ss)
    
    #with multiprocessing.Pool(processes=num_parallel) as pool:
    #    set_rings_all = set() # set to store enumerated rings
    #    set_rings = pool.map(fun_enum_rings, list_nodes)
    #    for s in set_rings:
    #        if len(s)>0:
    #            for ss in list(s):
    #                set_rings_all.add(ss)
    
    return set_rings_all

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
    progress_bar = tqdm(total = num)
    
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
        progress_bar.update()
        
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

def sub_parallel_enumerate_guttman_ring(i, atoms_extracted, G, chemical_bond_index):
    #print(f'process {i} done!')
    import copy
    
    _G = copy.deepcopy(G)
    
    set_rings = set()
        
    n0 = chemical_bond_index[i,0]
    n1 = chemical_bond_index[i,1]
    if (n0 in atoms_extracted)|(n1 in atoms_extracted):            
        if has_igraph:
            ei = _G.get_eid(n0, n1)
            _G.delete_edges(ei)
        else:
            _G.remove_edge(n0, n1)
        try:
            if has_igraph:
                paths = _G.get_all_shortest_paths(n0, n1)
            else:
                paths = nx.all_shortest_paths(_G, source=n0, target=n1)                    
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
            _G.add_edges([(n0, n1)])
        else:
            _G.add_edge(n0, n1)
        
    return set_rings

def parallel_enumerate_guttman_ring(atoms_extracted, atoms_all, chemical_bond_index, num_parallel=-1):

    # the number of processes in parallel computation
    if num_parallel > os.cpu_count():
        num_parallel = os.cpu_count()-2
    elif num_parallel == -1:
        num_parallel = os.cpu_count()-2
    
    index_bonds = np.arange(chemical_bond_index.shape[0])
    #list_index_bonds = np.array_split(index_bonds, num_parallel)    
    
    # initialize graph object
    if has_igraph:
        G = ig.Graph()
        G.add_vertices(len(atoms_all))
        G.add_edges(chemical_bond_index)
    else:
        G = nx.Graph()
        G.add_nodes_from(atoms_all)
        G.add_edges_from(chemical_bond_index)
            
    # make function with fixed parameters
    fun_enum_rings = partial(sub_parallel_enumerate_guttman_ring, atoms_extracted=atoms_extracted, G=G, chemical_bond_index=chemical_bond_index) 
    num = chemical_bond_index.shape[0]
    
    with tqdm(total=num) as progress_bar:
        
        with ThreadPoolExecutor(max_workers=num_parallel) as executor:
            #progress_bar = tqdm(total = num)
            set_rings_all = set()
            procs = []
            for n in index_bonds:
                proc = executor.submit(fun_enum_rings, n)            
                procs.append(proc)
                
            for f in futures.as_completed(procs):
                progress_bar.update(1)
                                
            executor.shutdown(wait=True)
            
            for p in procs:
                s= p.result()            
                if len(s)>0:
                    for ss in list(s):
                        set_rings_all.add(ss)
    
    return set_rings_all

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
    progress_bar = tqdm(total = num)
    
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
        progress_bar.update()
        
    return set_rings
    
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

    def __init__(self):
        pass

class Ring(object):
            
    def __init__(self, atoms, indexes=None):
        self.atoms = atoms
        self.indexes = indexes
        self._roundness = None
        self._roughness = None
        self.ellipsoid_lengths = None
    
    def __getitem__(self, i):
        return self.indexes[i]
    
    def __str__(self):
        return str(self.indexes)
    
    @property
    def close(self):
        bond_index = self.atoms.bonds.tolist()
        sum_shift = np.zeros(3, dtype='int8')
        lpath = list(self.indexes)
        for ns, ne in zip(lpath, lpath[1:]+lpath[:1]):
            i = bond_index[ns].index(ne)
            shift = self.atoms.shifts[ns][i]            
            sum_shift += shift    
        if not np.all(sum_shift == 0):
            return False
        else:
            return True
        
    @property
    def number(self):
        return len(self.indexes)
    
    def calc_svd(self):
        centres = self.atoms.norm_positions        
        i = self.indexes[0]
        xyz = [[0., 0., 0.]]
        for j in self.indexes[1:]:
            x = centres[j][0]-centres[i][0]+3.
            y = centres[j][1]-centres[i][1]+3.
            z = centres[j][2]-centres[i][2]+3.
            x = 2.*(x/2.-int(x/2.))-1.
            y = 2.*(y/2.-int(y/2.))-1.
            z = 2.*(z/2.-int(z/2.))-1.
            x, y, z = np.dot(self.atoms.volume.vectors, [x,y,z])
            xyz.append([x, y, z])        
        xyz = np.array(xyz)
        ring_centers = xyz.mean(axis=0)
        xyz0 = xyz - ring_centers
        svd_u, svd_s, svd_uh = np.linalg.svd(xyz0)
        self.ellipsoid_lengths = svd_s
    
    @property
    def roundness(self):
        if self.ellipsoid_lengths is None:
            self.calc_svd()
        if self._roundness is None:
            self._roundness = self.ellipsoid_lengths[1]/self.ellipsoid_lengths[0]        
        return self._roundness
        
    @property
    def roughness(self):
        if self.ellipsoid_lengths is None:
            self.calc_svd()
        if self._roughness is None:            
            self._roughness = self.ellipsoid_lengths[2]/np.sqrt(self.ellipsoid_lengths[0]*self.ellipsoid_lengths[1])        
        return self._roughness

    @property
    def over_boundary(self):
        bond_index = self.atoms.bonds.tolist()
        lpath = list(self.indexes)
        for ns, ne in zip(lpath, lpath[1:]+lpath[:1]):
            i = bond_index[ns].index(ne)
            shift = self.atoms.shifts[ns][i]
            if not np.all(shift == 0):
                return True
        return False
    
class RINGs(Molecule):
    
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
                
    atom_symbols : array (str) of shape (num_atoms,)
        List of atomic symbols
    
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
    
    def __init__(self, atoms):
        super().__init__()
        
        self.atoms = atoms
        self.num_atoms = self.atoms.number
        self.atom_symbols = np.array(self.atoms.elements)
        self.rings = []
            
    def calculate(self, ring_type, pair_atom_symbols, p_pair=0.3,
                  cutoff_size=24, num_parallel=0, atoms_extracted=None, chain=False):
                
        index_atoms = []
        #print(self.atoms.bonds)
        for index, neighs in enumerate(self.atoms.bonds):
            for nei in neighs:
                if index < nei:
                    index_atoms.append([index, nei])
        index_atoms = np.array(index_atoms)
                
        # judge if the distance of each pairs of atoms is less than pair_dist
        bool_bond = np.zeros(index_atoms.shape[0], dtype=bool)
        for (a0, a1) in (pair_atom_symbols):
            bool_bond = (((self.atom_symbols[index_atoms[:,0]]==a0) & (self.atom_symbols[index_atoms[:,1]]==a1)) \
                | ((self.atom_symbols[index_atoms[:,0]]==a1) & (self.atom_symbols[index_atoms[:,1]]==a0))) | bool_bond
        self.chemical_bond_index_atoms = index_atoms[bool_bond, :]
                
        self.enumerate_ring(ring_type, cutoff_size, num_parallel, atoms_extracted)
        
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
                
        flag_type = True
        if (ring_type == RINGs.RingType.PRIMITIVE)&(num_parallel == 0):            
            set_rings = enumerate_primitive_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, cutoff_size)
        elif (ring_type == RINGs.RingType.PRIMITIVE)&(num_parallel != 0):
            set_rings = parallel_enumerate_primitive_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, 
                                                          cutoff_size, num_parallel=num_parallel)
        elif (ring_type == RINGs.RingType.PRIMITIVE_KING)&(num_parallel == 0):
            set_rings = enumerate_king_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, flag_primitive=True)
        elif (ring_type == RINGs.RingType.KING)&(num_parallel == 0):
            set_rings = enumerate_king_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, flag_primitive=False)
        elif (ring_type == RINGs.RingType.KING)&(num_parallel != 0):
            set_rings = parallel_enumerate_king_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, \
                 flag_primitive=False, num_parallel=num_parallel)
        elif (ring_type == RINGs.RingType.GUTTMAN)&(num_parallel == 0):
            set_rings = enumerate_guttman_ring(atoms_extracted, atoms_all,self.chemical_bond_index_atoms)
        elif (ring_type == RINGs.RingType.GUTTMAN)&(num_parallel != 0):
            set_rings = parallel_enumerate_guttman_ring(atoms_extracted, atoms_all, self.chemical_bond_index_atoms, 
                                                        num_parallel=num_parallel)
        else:
            print('The indicated ring type is incorrect.')
            print('Choose a fring tyep from Primitive, King, Primitive_King, and Guttman.') 
            flag_type=False
        
        if flag_type:
            self.ring_type = ring_type
            self.rings = [Ring(self.atoms, list(r)) for r in list(set_rings)] # save rings as list
        