import numpy as np
from sidechain_cap_module import *

#SCRIPT SETTINGS
#-----------------------------------------------------------------------------------------------------------------------------------------------------

#control the cutoff below which you consider an atom pair chemically bonded
radial_cutoff = 1.80

#planarity used to determine aromatic atoms, choose the percentage buffer around the plane-bond angle needed for planarity (90 deg.)
#larger planarity tolerance --> angles with more deviation from 90 deg. will be considered planar
planarity_tolerance = 0.15

file_1 = 'mol1.xyz' #name of 1st input pos file
file_2 = 'mol3.xyz' #name of 2nd input pos file

#read in atomic coords and define atomic numbers
coords_1 = xyz_parser(file_1)
coords_2 = xyz_parser(file_2)

atom_key = {'S':16, 'N':7, 'C':6, 'O':8, 'H':1, 'F':9}
reverse_atom_key = {16:'S', 7:'N', 6:'C', 8:'O', 1:'H', 9:'F'}
#------------------------------------------------------------------------------------------------------------------------------------------

#re-writing nuclear coordinates with atomic numbers instead of element labels
coords_1 = replace_numbers(coords_1, atom_key)
coords_2 = replace_numbers(coords_2, atom_key)

mol1_nn, mol1_planarity = characterise_molecule(coords_1, radial_cutoff)
mol2_nn, mol2_planarity = characterise_molecule(coords_2, radial_cutoff)

mol1_sorted_nn = sort_nn_by_distance(coords_1, mol1_nn)
mol2_sorted_nn = sort_nn_by_distance(coords_2, mol2_nn)

rb_1, cb_1, tb_1 = get_bond_pairs(coords_1, mol1_sorted_nn, mol1_planarity)
rb_2, cb_2, tb_2 = get_bond_pairs(coords_2, mol2_sorted_nn, mol2_planarity)

all_bonds_1 = rb_1 + cb_1 + tb_1
all_bonds_2 = rb_2 + cb_2 + tb_2

mapping_dict = map_atom_indices(all_bonds_1, all_bonds_2, 184)

#re-writing the .xyz file of mol2 whose atomic indices are being rearranged to match 1
new_posfile = open(f'{file_2[:-4]}-mapped.xyz','a')

new_posfile.write(f'{len(coords_2)}' + '\n')
new_posfile.write(' ' + '\n')

for index in range(len(coords_1)):

    mol1_to_mol2_index = mapping_dict[index]

    row = coords_2[mol1_to_mol2_index]
    atomic_number = row[0]
    element = reverse_atom_key[atomic_number]

    new_posfile.write(f'{element}   ' + str(row[1:])[1:-1] + '\n')

new_posfile.close()
