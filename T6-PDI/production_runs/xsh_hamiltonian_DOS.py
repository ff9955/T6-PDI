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

    return eigenvals, eigenvecs

#CSS_indices = np.array([41, 42, 43, 44, 45, 46, 47, 48, 49, 90, 92, 93, 94, 95, 96, 97, 98, 99, 140, 141, 143, 144, 145, 146, 147, 148, 149, 190, 191, 192, 194, 195, 196, 197, 198, 199, 240, 241, 242, 243, 245, 246, 247, 248, 249, 346, 347, 348, 349, 445, 447, 448, 449, 545, 546, 548, 549, 645, 646, 647, 649, 745, 746, 747, 748])
CSS_indices = np.array([9])
#indices of CT-diabats beyond coulomb barrier peak

number_excitons=10
number_diabats=60

CSS_list = []
XT_list = []
CT_list = []
iCT_list = []
H_list = []

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
                    eigenvals, eigenvecs = build_diagonalise(parsed_H, number_diabats)

                    active_eigenvector = eigenvecs[:, active_state]

#                    if counter < 100:
#                        counter += 1
#                        parsed_H = []
#                        H_file.readline()
#                        continue
#                    elif counter > 100:
#                        break
#                    elif np.sum(active_eigenvector[-number_excitons:]**2) < 0.95:
#                        break

                    for index in range(len(eigenvecs)):

                        single_vector = eigenvecs[:,index]
                        single_energy = eigenvals[index]

                        if np.sum(single_vector[np.array([30])]**2) > 0.80:
                            iCT_list.append(single_energy)
                        #if single_vector[360]**2 > 0.80:
                        #   iCT_list.append(single_energy)
                        elif np.sum(single_vector[CSS_indices]**2) > 0.80:
                            CSS_list.append(single_energy)
                        elif np.sum(single_vector[:-number_excitons]**2) > 0.95:
                            CT_list.append(single_energy)
                        elif np.sum(single_vector[-number_excitons:]**2) > 0.95:
                            XT_list.append(single_energy)
                        else:
                            H_list.append(single_energy)
    
                    parsed_H = []
                    H_file.readline()
                    counter += 1

                    if counter > 0:
                        break #breaking here to sample DOS of first timestep

text_file_list = ['physopt_shortchain_Eb190_CT_nDOS.txt', 'physopt_shortchain_Eb190_XT_nDOS.txt', 'physopt_shortchain_Eb190_H_nDOS.txt','physopt_shortchain_Eb190_iCT_nDOS.txt', 'physopt_shortchain_Eb190_CSS_nDOS.txt']

energy_lists = [CT_list, XT_list, H_list, iCT_list, CSS_list]

for j in range(len(text_file_list)):

    data_file = open(text_file_list[j], 'a')
    for element in energy_lists[j]:
        data_file.write(str(element) + '\n')
    data_file.close()
