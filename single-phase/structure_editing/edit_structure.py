import numpy as np
from sidechain_cap_module import *

#This is a script used to cut aliphatic side chains off a core of aromatic rings. Upon passing the position file of a molecule
# (or multiple molecules), the script will read the coordinates of each molecule and notify the user when it has identified a side chain.
#It will then ask the user whether they want to cut the side chain and replace it with a single atom of the user's choice. The script currently 
#does not cap the side chains with anything other than single atoms. The input file of the molecule(s) must be in .xyz format. After cutting 
#the selected side chains, the script will save the new coordinates in a separate .xyz file.

#The side chains are identified by first distuingishing them from atoms inside aromatic rings. This is done by gathering the bonded atoms of each
#atom in the molecule and then measuring the planarity of each atom's environment. Without editing any of the source code, the script is controlled 
#by the set of parameters highlighted below. These will have to be changed for each molecule type for the script to run effectively.

#SCRIPT SETTINGS
#-----------------------------------------------------------------------------------------------------------------------------------------------------

#number of atoms per molecule
atoms_per_mol = 187
#number of molecules in the position file
number_molecules = 1

#control the cutoff below which you consider an atom pair chemically bonded
radial_cutoff = 1.85

#planarity used to determine aromatic atoms, choose the percentage buffer around the plane-bond angle needed for planarity (90 deg.)
#larger planarity tolerance --> angles with more deviation from 90 deg. will be considered planar
planarity_tolerance = 0.15

#choose the atom you want to replace the side chains
cap_atom = 'H'

file_name = 'mol1.xyz' #name of input pos file
tag = 'mol1' #label for output files

#read in atomic coords and define atomic numbers
all_coordinates = xyz_parser(file_name)
atom_key = {'S':16, 'N':7, 'C':6, 'O':8, 'H':1, 'F':9}
reverse_atom_key = {16:'S', 7:'N', 6:'C', 8:'O', 1:'H', 9:'F'}

#-----------------------------------------------------------------------------------------------------------------------------------------------------

#re-writing nuclear coordinates with atomic numbers instead of element labels
all_coordinates = replace_numbers(all_coordinates, atom_key)

#empty lists for storing bond lengths and angles in the order: aromatic C/N/S then aliphatic
total_bond_lists = [[], [], [], []]

#lists for storing internal bond angles in same order as above
internal_angle_lists = [[], [], [], []]

#iterating through all molecules in pos file containing only 1 molecule type
mol_counter = 0
for index in range(0, len(all_coordinates), atoms_per_mol):
    mol_counter += 1
    print(f'Molecule {mol_counter}')

    #indexing coordinates of single molecule in pos file
    coordinates = all_coordinates[index:index+atoms_per_mol]

    #getting list of each atom's nearest neighbours (below cutoff) and angles between bond vectors and their cross products
    total_nn, total_planarity = characterise_molecule(coordinates, radial_cutoff)

    #diving atoms into three categories depending on if they're in a ring, aliphatic chain, or terminal
    ring_bonds, chain_bonds, terminal_bonds = get_bond_pairs(coordinates, total_nn, total_planarity)

    #extracting indices of atoms in aromatic rings (i.e in a ring and a planar environment)
    aromatic_atom_indices, atom_type_list = get_aromatics(coordinates, total_nn, total_planarity, ring_bonds, planarity_tolerance)
    print(aromatic_atom_indices)

    #identify the length of each side chain and ask user whether it will be cut or not and replaced with the cap atom
    chains_to_cut = select_side_chains(aromatic_atom_indices, chain_bonds, terminal_bonds, cap_atom)

    #remove selected side chains and replace with cap atom
    new_coordinates = cut_chains(coordinates, total_nn, chains_to_cut, atom_key, cap_atom)

    #re-writing new file
    new_posfile = open(f'{tag}-{mol_counter}-capped.xyz','a')

    new_posfile.write(f'{len(new_coordinates)}' + '\n')
    new_posfile.write(' ' + '\n')

    for row in new_coordinates:
    
        atomic_number = row[0]
        element = reverse_atom_key[atomic_number]
    
        new_posfile.write(f'{element}   ' + str(row[1:])[1:-1] + '\n')

    new_posfile.close()

    #get bond lengths and internal bond angles of different atom types
    bl_list, ba_list = collect_bonds_angles(coordinates, atom_type_list, total_nn)

    for index in range(len(bl_list)):

        total_bond_lists[index].append(bl_list[index])
        internal_angle_lists[index].append(ba_list[index])

#aromatic C/N/S atoms, then aliphatic non-terminal atoms
#array_labels = ['arC', 'arN', 'arS', 'al']

#for index in range(len(total_bond_lists)):

#    np.savetxt(f'{tag}-{array_labels[index]}-bonds.txt', np.array(total_bond_lists[index]))
#    np.savetxt(f'{tag}-{array_labels[index]}-internal-angles.txt', np.array(internal_angle_lists[index]))
