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

    hamiltonian = hamiltonian + hamiltonian.T

    eigenvals, eigenvecs = np.linalg.eigh(hamiltonian)
    idx = eigenvals.argsort()

    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]

    return eigenvals, eigenvecs, hamiltonian

number_excitons=20
number_diabats=420

XT_relaxation_list = []
CT_relaxation_list = []
iCT_relaxation_list = []

iXT_site_energies = []
iCT_site_energies = []

for trajectory in trajectory_list:

    coulomb_barrier = np.loadtxt(f'{target_directory}/FRZ_DIAGONALS_COULOMB.include')
    coulomb_barrier = coulomb_barrier[:,-1]

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
                    eigenvals, eigenvecs, H = build_diagonalise(parsed_H, number_diabats)

                    nace_file.readline()
                    nace_file.readline()
                    parsed_nace = nace_file.readline().split()
                    parsed_nace = np.array(parsed_nace, dtype=float)

                    active_state = int(parsed_nace[0]-1)
 
                    site_energies = np.diag(H)
                    CT_energies = site_energies[:-number_excitons] - coulomb_barrier

                    active_eigenvec = eigenvecs[:,active_state]
                    CT_populations = active_eigenvec[:-number_excitons]**2
                    XT_populations = active_eigenvec[-number_excitons:]**2

                    #conditional for selecting site E of interfacial exciton, when IPR ~ 1
                    if np.sum(XT_populations) > 0.95:
                        if XT_populations[0] > 0.9:
                            iXT_site_energies.append(H[400,400])

                    #selecting site E if interfacial CT-state, when e/h IPR ~ 1
                    elif np.sum(XT_populations) < 0.05:
                        if CT_populations[360] > 0.9:
                            iCT_site_energies.append(H[360,360])

                    #XT_energies = site_energies[-number_excitons:] - CT_energies[360] 
                    XT_energies = site_energies[-number_excitons:]

                    CT_populations = CT_populations/np.sum(CT_populations)
                    XT_populations = XT_populations/np.sum(XT_populations)

                    exp_XT_energy = np.sum(XT_populations*XT_energies) - np.max(XT_energies)
                    XT_relaxation_list.append(exp_XT_energy)

                    exp_CT_energy = np.sum(CT_populations*CT_energies) - np.max(CT_energies)
                    CT_relaxation_list.append(exp_CT_energy)

                    iCT_energy = CT_energies[360] - np.max(CT_energies) 
                    iCT_relaxation_list.append(iCT_energy)

                    parsed_H = []
                    H_file.readline()
                    counter += 1

filename='physopt_retry'

#data_file = open(f'{target_directory}/{filename}_CT_site_energies.txt', 'a')
#data_file.write(trajectory_list[0] + ', ' + str(CT_relaxation_list)[1:-1] + '\n')
#data_file.close()

data_file = open(f'{target_directory}/{filename}_XT_site_energies.txt', 'a')
data_file.write(trajectory_list[0] + ', ' + str(XT_relaxation_list)[1:-1] + '\n')
data_file.close()

#data_file = open(f'{target_directory}/{filename}_iCT_site_energies.txt', 'a')
#data_file.write(trajectory_list[0] + ', ' + str(iCT_relaxation_list)[1:-1] + '\n')
#data_file.close()

#data_file = open(f'{filename}_localised_iXT_energies.txt', 'a')
#for element in iXT_site_energies:
#    data_file.write(str(element) + '\n')
#data_file.close()
#
#data_file = open(f'{filename}_localised_iCT_energies.txt', 'a')
#for element in iCT_site_energies:
#    data_file.write(str(element) + '\n')
#data_file.close()
