import sys
import math
import numpy as np
import numpy.linalg as la
from . import file
import collections
import itertools

def clamping_acos(cos):
    """
    Calculate arccos with its argument clamped to [-1, 1]
    """
    if cos > 1:
        return 0
    if cos < -1:
        return math.pi/2
    return math.acos(cos)

def get_bonds_with_constant_delta(atoms, delta):
    """
    Two atoms are connected if their distance is less than or equals a given delta.
    """
    bond_target_index_arrays = []
    for index, position in enumerate(atoms.positions):
        distance_vector_array = atoms.positions[index+1:] - position
        distance_squared_array = np.sum(np.square(distance_vector_array), axis=1)
        delta_squared = delta**2
        bond_target_indices = (distance_squared_array <= delta_squared).nonzero()[0] + index + 1
        bond_target_index_arrays.append(bond_target_indices)
    return bond_target_index_arrays

class MaximumRadiusException(Exception):
    pass

def get_bonds_with_radii(atoms, radii_sum_factor):
    """
    Two atoms are connected if their distance is less than or equals the sum of
    their covalent radii times a radii_sum_factor (e.g. 1.15).
    """
    bond_target_index_arrays = []
    length_arrays = []
    shift_arrays = []
    if not atoms.volume.periodic:
        for index, position in enumerate(atoms.positions):
            distance_vector_array = atoms.positions[index+1:] - position
            distance_squared_array = np.sum(np.square(distance_vector_array), axis=1)
            delta_squared = np.square((atoms.covalence_radii[index+1:] + atoms.covalence_radii[index]) * radii_sum_factor)            
            bond_target_indices = (distance_squared_array <= delta_squared).nonzero()[0] + index + 1
            bond_target_index_arrays.append(bond_target_indices)
            length_arrays.append(np.sqrt(distance_squared_array[bond_target_indices - index - 1]))       
            shift_arrays.append(np.array([[0, 0, 0] for _ in bond_target_indices]))
    else:
        dist_max = 0.0
        #for ri, rj in itertools.combinations(atoms.covalence_radii, 2):
        #    dist_max = max(dist_max, ri+rj)
        for length in atoms.bond_lengths.values():
            dist_max = max(dist_max, length)
        grid = atoms.grid
        
        # TODO : Considering whether to expand the cell
        try:        
            if dist_max > grid.rmax:
                error = "*** Error : Maximum distance exceeds cell size in function bonds.get_bonds_with_radii."
                raise MaximumRadiusException(error)
                #dist_max = grid.rmax
        except MaximumRadiusException as e:
            print(e)
            print("Use a supercell structure to ensure that the minimum cell size exceeds {:.2f}.".format(grid.rmax))
            return
        
        bond_matrix = np.identity(len(atoms.symbol_set))
        for i, elem1 in enumerate(atoms.symbol_set):
            for j, elem2 in enumerate(atoms.symbol_set):
                pair = (elem1, elem2)
                if pair in atoms.bond_lengths.keys():
                    bond_matrix[i][j] = atoms.bond_lengths[pair]
                else:
                    pair = (elem2, elem1)
                    bond_matrix[i][j] = atoms.bond_lengths[pair]
        elements_indices = [atoms.symbol_set.index(e) for e in atoms.symbols]
        for index in range(atoms.num_total):
            grid.neighbours(index, dist_max, index+1, atoms.num_total)
            #delta = (atoms.covalence_radii[grid.inei] + atoms.covalence_radii[index]) * radii_sum_factor
            delta = np.array([bond_matrix[elements_indices[index]][elements_indices[j]] for j in grid.inei]) * radii_sum_factor           
            bond_target_indices = np.array(grid.inei)
            indices = (grid.d <= delta).nonzero()[0]
            bond_target_index_arrays.append(bond_target_indices[indices])
            length_arrays.append(np.array(grid.d)[indices])
            shift_arrays.append(np.array(grid.shifts).T[indices])

    # get_bonds_symetric_indicies
    bond_target_index_arrays = [bond_target_index_array.tolist() for bond_target_index_array in bond_target_index_arrays]
    length_arrays = [length_array.tolist() for length_array in length_arrays]
    shift_arrays = [shift_array.tolist() for shift_array in shift_arrays]
    for index, bond_target_index_array in enumerate(bond_target_index_arrays):
        for i, bond_target_index in enumerate(bond_target_index_array):
            if bond_target_index > index:
                bond_target_index_arrays[bond_target_index].append(index)
                length_arrays[bond_target_index].append(length_arrays[index][i])
                shift_arrays[bond_target_index].append(-np.array(shift_arrays[index][i]))                

    bond_target_index_arrays = np.array(bond_target_index_arrays, dtype=object)  
    atoms.lengths = np.array(length_arrays, dtype=object)
    atoms.shifts = [] # np.array(shift_arrays, dtype=object)
    for shift in shift_arrays:
        atoms.shifts.append(np.array(shift))
        
    return bond_target_index_arrays

def get_bonds_with_element_pair(atoms, min_dist, max_dist, radii_sum_factor):
    bond_target_index_arrays = []
    for index, position in enumerate(atoms.positions):
        distance_vector_array = atoms.positions[index+1:] - position
        distance_squared_array = np.sum(np.square(distance_vector_array), axis=1)
        ntypes = len(atoms.symbol_set)
        #ic = int(i*(2*ntypes-(i+1))/2+j)
        ics = atoms.indices[index+1:]*(2*ntypes-(atoms.indices[index+1:]+1))/2+atoms.indices[index]
        ics.astype(int)
        bonds = np.zeros_like(ics,dtype=float)
        for ic, dist in max_dist.items():
            bonds[ics == ic] = dist
        delta_squared = np.square(bonds * radii_sum_factor)
        bond_target_indices = (distance_squared_array <= delta_squared).nonzero()[0] + index + 1
        bond_target_index_arrays.append(bond_target_indices.tolist())
    return bond_target_index_arrays
        
def get_bonds_symetric_indicies(bond_target_index_arrays):
    """
    inserts all symetric bonds into the of bonds for each atom
    :param bond_target_index_arrays: list with all symetric bond indices
    :return:
    """
    bond_target_index_arrays = [bond_target_index_array.tolist() for bond_target_index_array in bond_target_index_arrays]
    for index, bond_target_index_array in enumerate(bond_target_index_arrays):
        for bond_target_index in bond_target_index_array:
            if bond_target_index > index:
                bond_target_index_arrays[bond_target_index].append(index)

    bond_target_index_arrays = np.array(bond_target_index_arrays, dtype=object)    
    return bond_target_index_arrays

def calculate_bond_length(atoms, i, j):
    pass

def calculate_bond_angles(atoms, bond_target_index_arrays):
    # Build dict: atom index -> bonded atom indices
    atom_bonds = collections.defaultdict(list)
    for source_index, target_indices in enumerate(bond_target_index_arrays):
        for target_index in target_indices:
            atom_bonds[source_index].append(target_index)
            atom_bonds[target_index].append(source_index)

    # Calculate angles between bonds sharing an atom
    bond_angles = {}
    for shared_atom_index, target_indices in atom_bonds.items():
        shared_atom_position = atoms.positions[shared_atom_index]
        for i, target_index1 in enumerate(target_indices):
            vec1 = atoms.positions[target_index1]-shared_atom_position
            nvec1 = vec1/np.linalg.norm(vec1)
            for target_index2 in target_indices[i+1:]:
                vec2 = atoms.positions[target_index2]-shared_atom_position
                nvec2 = vec2/np.linalg.norm(vec2)
                bond_angle = clamping_acos(np.dot(nvec1, nvec2))
                bond_angles[((target_index1, shared_atom_index), (shared_atom_index, target_index2))] = bond_angle
                bond_angles[((target_index2, shared_atom_index), (shared_atom_index, target_index1))] = bond_angle

    # Calculate bond chains of length 3
    bond_chains = []
    for source_index in atom_bonds.keys():
        bond_chains += find_bond_chains(atoms, source_index, atom_bonds)
    bond_chains_without_duplicates = set()
    for bond_chain in bond_chains:
        if tuple(reversed(bond_chain)) not in bond_chains_without_duplicates:
            bond_chains_without_duplicates.add(tuple(bond_chain))
    bond_chains = bond_chains_without_duplicates

    # Calculate angle in these chains around the axis of the connecting bond
    # according to http://en.wikipedia.org/wiki/Dihedral_angle
    bond_chain_angles = {}
    for bond_chain in bond_chains:
        index1, index2, index3, index4 = bond_chain
        axis = normalized(atoms.positions[index3] - atoms.positions[index2])

        vec1 = atoms.positions[index1] - atoms.positions[index2]
        vec1 -= np.dot(vec1, axis)*axis
        nvec1 = normalized(vec1)

        vec2 = atoms.positions[index4] - atoms.positions[index3]
        vec2 -= np.dot(vec2, axis)*axis
        nvec2 = normalized(vec2)

        angle = clamping_acos(np.dot(nvec1, nvec2))
        if np.dot(nvec2, np.cross(nvec1, axis)) < 0:
            angle = -angle
        bond_chain_angles[tuple(bond_chain)] = angle

    #print "bond_anglges", bond_angles[0]
    return bond_angles, bond_chain_angles

def normalized(vec):
    """ Return the normalized version of a numpy array """
    return vec/la.norm(vec)


def find_bond_chains(atoms, source_index, atom_bonds, length=3, previous_bond_chain=None):
    if previous_bond_chain is None:
        previous_bond_chain = [source_index]
    length -= 1
    bond_chains = []
    for target_index in atom_bonds[source_index]:
        if target_index not in previous_bond_chain:
            new_bond_chain = previous_bond_chain + [target_index]
            new_source_index = target_index
            if length > 0:
                new_bond_chains = find_bond_chains(atoms, new_source_index, atom_bonds, length, new_bond_chain)
                if new_bond_chains:
                    bond_chains += new_bond_chains
            else:
                bond_chains.append(new_bond_chain)
    return bond_chains

def export_bonds(filename, atoms):
    bond_target_index_arrays = get_bonds_with_constant_delta(atoms, 2.8)

    with open(filename, "w") as outfile:
        for source_index, target_indices in enumerate(bond_target_index_arrays):
            for target_index in target_indices:
                outfile.write("{} {}\n".format(source_index, target_index))

def export_bond_angles(filename, atoms):
    bond_target_index_arrays = get_bonds_with_constant_delta(atoms, 2.8)
    bond_angles, bond_chain_angles = calculate_bond_angles(atoms, bond_target_index_arrays)

    with open(filename, "w") as outfile:
        for bond1, bond2 in bond_angles.keys():
            if bond1[0] > bond2[1]:
                outfile.write("{} {} {} {}\n".format(bond1[0], bond1[1], bond2[1], bond_angles[bond1, bond2]))

def export_bond_dihedral_angles(filename, atoms):
    bond_target_index_arrays = get_bonds_with_constant_delta(atoms, 2.8)
    bond_angles, bond_chain_angles = calculate_bond_angles(atoms, bond_target_index_arrays)

    with open(filename, "w") as outfile:
        for bond_chain, angle in bond_chain_angles.items():
            outfile.write("{} {} {} {}".format(*bond_chain))
            outfile.write(" {}\n".format(angle))

