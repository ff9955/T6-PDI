import numpy as np
from file_parsers import *
import sys
from numba import jit

target_directory = sys.argv[1]
trajectory_list = [sys.argv[2]]

@jit
def build_diagonalise(pseudoH_lines, number_diabats):

    hamiltonian = np.zeros((number_diabats, number_diabats))

    for element in pseudoH_lines:

       row = int(element[0]) - 1
       column = int(element[1]) - 1
       quantity = element[2]

       if row == column:
           quantity = quantity/2

       hamiltonian[row,column] = quantity

    hamiltonian = hamiltonian - np.identity(number_diabats)*hamiltonian[30,30]

    hamiltonian = hamiltonian + hamiltonian.T

    eigenvals, eigenvecs = np.linalg.eigh(hamiltonian)
    idx = eigenvals.argsort()

    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]

    return eigenvals, eigenvecs, hamiltonian

number_excitons=10
number_diabats=60

energy_list = []
band_rank_list = []

for trajectory in trajectory_list:

    with open(f'{target_directory}/run-fssh-{trajectory}/run-pseudo-hamilt-1.xyz') as H_file:
        with open(f'{target_directory}/run-fssh-{trajectory}/run-nace-active-1.xyz') as nace_file:

            H_file.readline()
            H_file.readline()
            parsed_H = []
            counter = 0

            for line in H_file:

                split_line = line.split()
                if len(split_line) == 3:
                    parsed_H.append(split_line)

                else:
                    parsed_H = np.array(parsed_H, dtype=float)

                    nace_file.readline()
                    nace_file.readline()
                    parsed_nace = nace_file.readline().split()
                    parsed_nace = np.array(parsed_nace, dtype=float)

                    active_state = int(parsed_nace[0]-1)
                    eigenvals, eigenvecs, H = build_diagonalise(parsed_H, number_diabats)
                    active_energy = eigenvals[active_state] #indexing active state energy

                    #ict_energy = H[360,360]

                    CT_list = []
                    XT_list = []
                    hybrid_list = [] #initialising lists of energies of different eigenstate types

                    for index2 in range(len(eigenvecs)):

                        single_vector = eigenvecs[:,index2]
                        single_energy = eigenvals[index2]

                        if np.sum(single_vector[:-number_excitons]**2) > 0.95:
                            CT_list.append(single_energy) #if state is CT, pass into CT energy list, etc...
                        elif np.sum(single_vector[-number_excitons:]**2) > 0.95:
                            XT_list.append(single_energy)
                        else:
                            hybrid_list.append(single_energy)

                    #I define adiabat energy wrt the lower edge of each band or the interfacial CT diabat
                    if active_energy in CT_list:
                        energy_rank = CT_list.index(active_energy) + 1
                        relative_energy = active_energy - np.min(CT_list)
#                        relative_energy = active_energy - ict_energy
                    elif active_energy in XT_list:
                        energy_rank = XT_list.index(active_energy) + 1
                        relative_energy = active_energy - np.min(XT_list)
#                        relative_energy = active_energy - ict_energy
                    else:
                        energy_rank = hybrid_list.index(active_energy) + 1
                        relative_energy = active_energy - np.min(hybrid_list)
#                        relative_energy = active_energy - ict_energy

                    band_rank_list.append(energy_rank)
                    energy_list.append(relative_energy)

                    parsed_H = []
                    H_file.readline()

filename='physopt_shortchain_Eb190'

data_file = open(f'{target_directory}/{filename}_active_Erank.txt', 'a')
data_file.write(trajectory_list[0] + ', ' + str(band_rank_list)[1:-1] + '\n')
data_file.close()

data_file = open(f'{target_directory}/{filename}_active_energy.txt', 'a')
data_file.write(trajectory_list[0] + ', ' + str(energy_list)[1:-1] + '\n')
data_file.close()


