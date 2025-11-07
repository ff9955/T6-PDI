import numpy as np
from scipy.spatial import KDTree, transform
import sys
#from collections import Counter

#python module containing functions to delete side chains of an aromatic core and replace them with single terminal atoms
#this is done by first gathering a group of nearest neighbour for each atom then dividing these into aromatic, aliphatic or terminal 
#atoms

def xyz_parser(file_path):
#simple parser to get coordinates of single frame
    
    with open(file_path) as full_file:

        all_lines = full_file.readlines()
        all_lines = all_lines[2:]
        
        split_lines = map(lambda x: x.split(), all_lines)
        split_lines = list(split_lines)

    output_array = np.array(split_lines)
    return output_array

def get_bond_angle(vector1, vector2):
    #get angle between 2 bonds using scalar product

    v1_norm = np.linalg.norm(vector1)
    v2_norm = np.linalg.norm(vector2)
    dp = np.dot(vector1, vector2)

    cos_theta = dp/(v1_norm*v2_norm)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    angle = np.arccos(cos_theta)
    angle = np.degrees(angle)

    return angle

def divide_atom_types(number_atoms, bond_pairs, total_planarity):
    #given an input of the indices of atoms that are directly bound, sort these pairs into bonds that are in rings, chains or terminal 
    #groups

    loop_done = False
    first_loop = True

    ring_atoms = []
    chain_atoms = []
    terminal_atoms = []

    #while loop, iterating over all atoms each loop
    while not loop_done:
        atom_freq = []

        #indexing each atom in molecule
        for atom_index in range(number_atoms):
            counter = 0

            relevant_pairs = []
            #indexing over all bond pairs
            for pair_index in range(len(bond_pairs)):

                pair = bond_pairs[pair_index]
                #if indexed atom is in the pair, add to counter and keep this pair in memory
                if atom_index in pair:
                    counter += 1
                    relevant_pairs.append(pair)
        
            #store the number of times the atom appeared in the bond list
            atom_freq.append(counter)

            #if you're in the first loop and the planarity has an actual value (i.e not H or terminal atom),
            #remove the bond pair from the total list of bond pairs and add atom to list of terminal atoms
            if (first_loop) and (not total_planarity[atom_index]):
                for pair in relevant_pairs:
                    bond_pairs.remove(pair)
                    terminal_atoms.append(pair)

            #this is only accessed after the first loop where the terminal atoms of the full molecule are removed, 
            #so a counter of 1 here means these are actually chain atoms
            elif counter == 1:
                for pair in relevant_pairs:
                    bond_pairs.remove(pair)
                    chain_atoms.append(pair)
            #we remove the bond pairs of terminal atoms from the full list each iteration in the loop

        first_loop = False
        #the loop keeps going until all terminal atoms have been deleted, leaving only ring atoms' bond pairs in the list
        if 1 not in atom_freq: loop_done = True

    #iterating through the remaining bond pairs list, which contains only ring atoms
    for pair in bond_pairs:
        #adding all atomic indices to ring atoms list
        ring_atoms.append(pair)
        #ring_atoms = ring_atoms + pair
    
    #removing duplicate ring atom indices, keeping the output as list
    #ring_atoms = list(set(ring_atoms))

    return ring_atoms, chain_atoms, terminal_atoms

def find_connection(query_atoms, bond_pairs):
    #given a list of indices of query atoms, find all atomic indices that are connected to them via bonds

    connected_atoms = []
    connection_found = False

    #for each query atom do the following
    for qa in query_atoms:
        print('atom:', qa)

        relevant_pairs = []
        for pair in bond_pairs:

            #isolate bond pairs with the query atom (qa) inside them
            if qa in pair:

                relevant_pairs.append(pair)
                #use the boolean to remember that you found a connection
                connection_found = True

                #remove qa from the relevant pairs, leaving only indices of atoms bound to qa
                pair.remove(qa)
                if pair[0] not in connected_atoms: connected_atoms.append(pair[0])

        print(relevant_pairs)
        
        #remove the relevant pairs from the total bond pair list
        for pair in relevant_pairs:
            bond_pairs.remove(pair)

    #if qa is connected to anything do the following
    if connection_found:
        #we have the indices to which qa is directly connected, so now find the indices to which the connected atoms are directly connected
        #since it's a recursive call, this will keep going until you find all possible connected atoms
        next_connected_atoms = find_connection(connected_atoms, bond_pairs)

        #add all connected atoms' indices to output list
        if next_connected_atoms:
            
            next_unique_connected_atoms = [element for element in next_connected_atoms if element not in connected_atoms]
            connected_atoms = connected_atoms + next_unique_connected_atoms

    #if qa has no connection, return False to end the recursion
    else:
        return False

    return connected_atoms

def replace_numbers(coordinates, atom_key):

    #read in atomic coords and define atomic numbers
    atom_key = {'S':16, 'N':7, 'C':6, 'O':8, 'H':1, 'F':9}

    #replace element symbol with atomic number
    for row in range(len(coordinates)):
        atom = coordinates[row][0]

        coordinates[row][0] = atom_key[atom]

    #converting coord array to array of floats
    vector_float = np.vectorize(float)
    coordinates = vector_float(coordinates)

    return coordinates

def characterise_molecule(coordinates, radial_cutoff):

    #getting each atom's nearest neighbour list within radial cutoff
    tree = KDTree(coordinates[:,1:])
    total_nn = tree.query_ball_point(coordinates[:,1:], r = radial_cutoff)
    #total_nn = list of lists of all atom's neighbours within the radial cutoff (defined above)

    #removing each instance of atom referencing itself in each neighbour list
    for index in range(len(total_nn)):
        total_nn[index].remove(index)

    #now onto identifying which atoms are aromatic
    total_planarity = []

    for index in range(len(total_nn)):

        #checking that the atom isn't hydrogen and that it's bound to >1 atom
        if (coordinates[index,0] != 1) and (len(total_nn[index]) > 1):
            nn = total_nn[index]

            #getting all bond vectors of given atom
            bond_vectors = [(coordinates[index,1:] - coordinates[neighbour,1:]) for neighbour in nn]
            #getting angles of all atom's bonds
            angles = [get_bond_angle(bond_vectors[0], bond_vectors[neighbour]) for neighbour in range(1,len(bond_vectors))]

            #getting vectors perpendicular to bond vectors
            cross_products = [np.cross(bond_vectors[0], bond_vectors[neighbour]) for neighbour in range(1,len(bond_vectors))]
            #testing planarity by getting angle between bond vectors and cross products of different bond vectors (should be close to 90 if planar)
            planarities = [get_bond_angle(cross_products[bond_pair], bond_vectors[-1 - bond_pair]) for bond_pair in range(len(cross_products))]

            total_planarity.append(np.mean(planarities))

        else:
            total_planarity.append(None)

    return total_nn, total_planarity
    
def get_bond_pairs(coordinates, total_nn, total_planarity):

    bond_pairs = []

    number_atoms = len(coordinates)

    #build list of indices of atoms that are chemically bonded
    for atom_index in range(number_atoms):

        #atom's neighbour list
        atoms_nn = total_nn[atom_index]

        for neighbour in atoms_nn:

            #append all pairs of neighbouring atoms if not already in list
            if atom_index < neighbour:
                bond_pairs.append([atom_index, neighbour])
            elif neighbour < atom_index:

                if [neighbour, atom_index] not in bond_pairs:
                    bond_pairs.append([atom_index, neighbour])

    #divide the bond pairs into those within rings, chains or at the end of chains
    ring_atoms, chain_atoms, terminal_atoms = divide_atom_types(number_atoms, bond_pairs, total_planarity)

    return ring_atoms, chain_atoms, terminal_atoms

def get_aromatics(coordinates, total_nn, total_planarity, ring_atoms, planarity_tolerance):

    all_aromatic_atoms = []
    aromatic_c_atoms = []
    aromatic_n_atoms = []
    aromatic_s_atoms = []
    aliphatic_atoms = []

    number_atoms = len(coordinates)

    #setting angle range within which we accept that molecule is planar (defined above)
    tol = 90*planarity_tolerance
    for atom_index in range(number_atoms):

        nn = total_nn[atom_index]
        planarity = total_planarity[atom_index]

    #1st check, must have planarity of aromatic atom (within tolerance), also excludes H or singly bound atoms, as they have planarity=None
        if not planarity:
            continue
        elif (planarity > (90 + tol)) or (planarity < (90 - tol)):
            aliphatic_atoms.append(atom_index)
            continue
        elif atom_index not in [element for sub in ring_atoms for element in sub]:
            aliphatic_atoms.append(atom_index)
            continue

    #checking that for aromatic atom, it's connected to 2 more atoms in planar environments (ring)
        planar_count = 0
        for neighbour in nn:
            planarity = total_planarity[neighbour]

            if not planarity:
                continue
            elif (planarity > (90 + tol)) or (planarity < (90 - tol)):
                continue

            planar_count += 1

        if (planar_count < 2):
            aliphatic_atoms.append(atom_index)
            continue

        all_aromatic_atoms.append(atom_index)

        if coordinates[atom_index,0] == 16:
            aromatic_s_atoms.append(atom_index)
        elif coordinates[atom_index,0] == 7:
            aromatic_n_atoms.append(atom_index)
        else:
            aromatic_c_atoms.append(atom_index)

    return all_aromatic_atoms, [aromatic_c_atoms, aromatic_n_atoms, aromatic_s_atoms, aliphatic_atoms]


def select_side_chains(all_aromatic_atoms, chain_atoms, terminal_atoms, cap_atom):

    chains_to_cut = []
    for atom_index in all_aromatic_atoms:

        #atom_index = index of aromatic atom
        #now in position to search for aliphatic side chains
        for pair in chain_atoms:

            #chain bond pair with aromatic atom in it --> side chain linked to aromatic core
            if atom_index in pair:

                #for this bond pair, see if aromatic atom is 1st or 2nd position
                oa_index = pair.index(atom_index)
                #use that info to get the index of the chain atom
                first_cut = pair[oa_index-1]

                #found only first atom of the side chain, now find all the others
                side_chain = find_connection([atom_index], chain_atoms)

                decide_cut = input(f'Found side chain with {len(side_chain)} atoms, cap with {cap_atom}? (y/n)')
                assert(decide_cut == 'y' or decide_cut == 'n'), 'Input should be y/n'

                if decide_cut == 'y':

                    #finding index of first side chain atom inside the list of all side chain atoms' indices
                    fc_index = side_chain.index(first_cut)
                    #make sure it's in the first index
                    if fc_index > 0: side_chain[0], side_chain[fc_index] = side_chain[fc_index], side_chain[0]
                    side_chain_terminals = []

                    #iterate through non-terminal indices of side chain atoms
                    for chain_atom in side_chain:
                        #if connected to terminal/H atom...
                        for pair in terminal_atoms:

                            if chain_atom in pair:

                                #add terminal/H atom index to list of atoms to cut
                                oa_index = pair.index(chain_atom)
                                side_chain_terminals.append(pair[oa_index-1])

                    #full chain of atoms to cut
                    side_chain = side_chain + side_chain_terminals

                    #append full chain to list of chains chosen to cut
                    chains_to_cut.append(side_chain)

    return chains_to_cut

def cut_chains(coordinates, total_nn, chains_to_cut, atom_key, cap_atom):

    number_atoms = len(coordinates)

    atoms_to_keep = []
    atoms_to_cap = []
    extra_cap_atoms = []

    for chain in chains_to_cut:
        #capping first atom of side chain with chosen atom (defined above as cap_atom)
        atoms_to_cap.append(chain[0])

    #keeping atoms that don't appear in the side chains you want to cut
    for atom_index in range(number_atoms):

        cut = False
        for chain in chains_to_cut:

            if atom_index in chain: cut = True
    
        if cut == False: atoms_to_keep.append(atom_index)

    #stopping here if no side chains selected
    if not atoms_to_cap:
        print('No chains to cap')
        return None
    

    #compiled all atoms to delete and which ones to replace
    #atoms_to_keep = np.array(atoms_to_keep)
    #atoms_to_cap = np.array(atoms_to_cap)

    #reassigning the atoms to cap with the atomic number of cap_atom or a methyl group
    if cap_atom=='Methyl':
        coordinates[atoms_to_cap,0] = atom_key['C'] #if methyl, assign 1st cap atom to C

        #go through all cap atoms indices
        for cap_index in atoms_to_cap:

            #get their neighbours
            nn = total_nn[cap_index]
            for n_index in nn:
                if n_index not in atoms_to_keep: #for neighbours of 1st cap atom that are not aromatic, replace with H
                    coordinates[n_index,0] = atom_key['H']
                    extra_cap_atoms.append(n_index)
    else:
        coordinates[atoms_to_cap,0] = atom_key[cap_atom]

    #deleting the rest from the molecule coordinates
    #atoms_to_keep = np.hstack((atoms_to_keep, atoms_to_cap))

    new_coordinates = coordinates[atoms_to_keep,:]
    if atoms_to_cap: new_coordinates = np.vstack((new_coordinates, coordinates[atoms_to_cap,:]))
    if extra_cap_atoms: new_coordinates = np.vstack((new_coordinates, coordinates[extra_cap_atoms,:]))

    return new_coordinates

def collect_bonds_angles(coordinates, atom_type_list, total_nn):
    #function to collect bond lengths and angles between all atoms in atom_type_list and their nearest neighbours

    #list of lists of bond lengths/angles for different elements, separate list in atom_type_list for each element
    bond_length_list = []
    bond_angle_list = []

    for atom_type in atom_type_list:

        #list for all bonds/angles involving this element
        bond_lengths = []
        bond_angles = []

        for atom_index in atom_type:

            nn = total_nn[atom_index]
            for n in nn:

                #neglecting bonds from c to heteratom as that will be sampled from the heteroatom lists
                if (atom_index in atom_type_list[0]) and (n not in atom_type_list[0]): continue

                #euclidean distance between neighbours
                bond_length = np.sqrt(np.sum((coordinates[atom_index,1:] - coordinates[n,1:])**2))
                bond_lengths.append(bond_length)

            #getting all bond vectors of given atom
            bond_vectors = [(coordinates[atom_index,1:] - coordinates[neighbour,1:]) for neighbour in nn]
            #getting angles of all atom's bonds

            for i in range(len(bond_vectors)-1):

                local_angles = [get_bond_angle(bond_vectors[i], bond_vectors[neighbour]) for neighbour in range(i+1,len(bond_vectors))]
                bond_angles = bond_angles + local_angles
            
        bond_length_list.append(bond_lengths)
        bond_angle_list.append(bond_angles)
    
    return bond_length_list, bond_angle_list

def sort_nn_by_distance(coordinates, total_nn):

    total_sorted_nn = []

    for atom_index in range(len(total_nn)):
        atom_coord = coordinates[atom_index, 1:]

        distances = []

        for nn in total_nn[atom_index]:
            
            neighbour_coord = coordinates[nn, 1:]
            distance = np.sqrt(np.sum((atom_coord - neighbour_coord)**2))
            distances.append(distance)
        
        distances = np.array(distances)
        idx = distances.argsort()
        sorted_nn = np.array(total_nn[atom_index])[idx]
        total_sorted_nn.append(sorted_nn)
    
    return total_sorted_nn

def map_atom_indices(mol_bonds_1, mol_bonds_2, common_atom):

    print('mol1 connections')
    mol_connections_1 = find_connection([common_atom], mol_bonds_1)
    print('mol3 connections')
    mol_connections_2 = find_connection([common_atom], mol_bonds_2)

    mapping_dict = dict()
    for index in range(len(mol_connections_1)):
        mapping_dict[mol_connections_1[index]] = mol_connections_2[index]

    mapping_dict[common_atom] = common_atom

    return mapping_dict
