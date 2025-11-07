import numpy as np
import sys

target_directory = sys.argv[1]

trajectory_range = [t for t in range(int(sys.argv[2]),int(sys.argv[3]))]

active_donor_molecules = 10
total_donor_molecules = 300
active_acceptor_molecules = 10
#specify number of electronically active donor/acceptor molecules, and total number of donors in the system

donor_atoms = 44
acceptor_atoms = 40
#specify number of atoms per single donor/acceptor molecule

atomic_weights = {'H':1.007825, 'O':15.994915, 'S':31.972072, 'N':14.000109, 'C':12}
#specify the atomic mass of all atoms of your donor/acceptor molecules are comprised, the script needs this to calculate the centres of mass of the molecules in each CT-state

total_donor_xCOMS = np.zeros((len(trajectory_range), active_donor_molecules))
total_acceptor_xCOMS = np.zeros((len(trajectory_range), active_acceptor_molecules))

with open(f'{target_directory}/DECOMP.include') as decomp_file:
    #DECOMP file contains the indices of the electronically active molecules; if this file is called something else, change the name here

    decomp_lines = decomp_file.readlines()
    for line in decomp_lines:

        split_line = line.split()
        if split_line[0] == 'INDEX_MOL_DECOMP':
            decomp_donor = list(map(int,split_line[1:active_donor_molecules + 1]))
            decomp_acceptor = list(map(int,split_line[active_donor_molecules + 1:]))
            decomp_indices = zip(decomp_donor, decomp_acceptor)
            decomp_enumerator = tuple(enumerate(decomp_indices))
            #decomp enumerator is a tuple of two iterables, the first being the indices of the elements in decomp_indices, and the second being the two-element lists themselves, which refer the indices of the active molecules in the position file
            break
       #parsing the DECOMP file to get separate lists of donor/acceptor indices


for trajectory in trajectory_range:
    #trajectory numbers used as indices, so the trajectories you read in must start at zero and be in ascending order, otherwise you will be indexing the wrong trajectories

    with open(f'{target_directory}/run-fssh-{trajectory}/pos-init.xyz') as position_file:

        position_lines = position_file.readlines()
        position_coords = list(map(lambda x: x.split(), position_lines))

        position_coords = position_coords[2:]
        position_coords = np.array(position_coords)
        #parsing the intial position file to get an array of every atom's xyz coordinates, and atom types

    atom_column = position_coords[:,0]
    #get 1D array of all atom types

    vector_map_fn = np.vectorize(lambda x: atomic_weights[x])
    atom_column = vector_map_fn(atom_column)
    #map strings of atom types onto their corresponding atomic masses

    coordinate_array = np.zeros(np.shape(position_coords))
    coordinate_array[:,0] = atom_column
    coordinate_array[:,1:] = position_coords[:,1:]
    #re-define coordinate array, but now with atomic masses instead of strings of atom types

    donor_COMS = np.zeros((active_donor_molecules, 3))
    acceptor_COMS = np.zeros((active_acceptor_molecules, 3))
    #initialise arrays for the COM of all active molecule, split into xyz components

    for index_number,index_pair in decomp_enumerator:
        #getting indices and DECOMP indices of donor/acceptor molecules

        donor_line = donor_atoms*(index_pair[0]-1)
        acceptor_line = donor_atoms*(total_donor_molecules) + acceptor_atoms*(index_pair[1] - 1 - total_donor_molecules)
        #getting the row corresponding to the first atom of the chosen donor molecule, and that of the acceptor molecule as well

        donor_coords = coordinate_array[donor_line: donor_line + donor_atoms]
        acceptor_coords = coordinate_array[acceptor_line: acceptor_line + acceptor_atoms]
        #use these initial lines and number_atoms per molecule to get the chunks of coordinates that belong to the chosen molecules

        total_donor_mass = np.sum(donor_coords[:,0])
        total_acceptor_mass = np.sum(acceptor_coords[:,0])
        #get molecular mass of each
    
        donor_masses = donor_coords[:,0][:,np.newaxis]
        donor_coords[:,1:] = np.multiply(donor_masses, donor_coords[:,1:])
        #multiply each atomic coordinates by atomic mass
        donor_com = np.sum(donor_coords[:,1:], axis=0)/total_donor_mass
        #sum weighted coordinates by different components, then divide by molecular mass to get COM, this gives you the xyz components of the COM

        acceptor_masses = acceptor_coords[:,0][:,np.newaxis]
        acceptor_coords[:,1:] = np.multiply(acceptor_masses, acceptor_coords[:,1:])

        acceptor_com = np.sum(acceptor_coords[:,1:], axis=0)/total_acceptor_mass
        #then just do the same for the acceptor

        donor_COMS[index_number] = donor_com
        acceptor_COMS[index_number] = acceptor_com
        #the chosen molecule pair's COMs are 3-element rows, which are then placed into arrays which will contains the xyz-COMS of all active donor and acceptor molecules

    total_donor_xCOMS[trajectory, :] = donor_COMS[:, 0]
    total_acceptor_xCOMS[trajectory, :] = acceptor_COMS[:, 0]

donor_xCOM_sd = np.std(total_donor_xCOMS, axis=0)
acceptor_xCOM_sd = np.std(total_acceptor_xCOMS, axis=0)

donor_xCOM_diffs = abs(np.diff(total_donor_xCOMS, axis=1))
avg_donor_xCOM_diffs = np.mean(donor_xCOM_diffs)

acceptor_xCOM_diffs = abs(np.diff(total_acceptor_xCOMS, axis=1))
avg_acceptor_xCOM_diffs = np.mean(acceptor_xCOM_diffs)

donor_xCOM_sd = donor_xCOM_sd/avg_donor_xCOM_diffs
acceptor_xCOM_sd = acceptor_xCOM_sd/avg_acceptor_xCOM_diffs

np.savetxt(f'{target_directory}/T6_initial_xCOMS.txt', total_donor_xCOMS)
np.savetxt(f'{target_directory}/PDI_initial_xCOMS.txt', total_acceptor_xCOMS)
