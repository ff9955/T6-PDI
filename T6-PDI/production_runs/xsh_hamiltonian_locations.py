import numpy as np
from file_parsers import xyz_parser, sh_log_parser
from xsh_analysis_functions import diabat_locations, IPR
import time
from numba import jit
import sys

target_directory = sys.argv[1]
trajectory = sys.argv[2]

number_diabats = 420
number_excitons = 20

trajectory_donor_xCOMS = np.loadtxt(f'{target_directory}/T6_initial_xCOMS.txt')
trajectory_acceptor_xCOMS = np.loadtxt(f'{target_directory}/PDI_initial_xCOMS.txt')

@jit
def inner_loop(pseudoH_lines, active_state_array, number_timesteps, number_diabats):

    file_length = len(pseudoH_lines)
    final_position_index = np.min(np.where(pseudoH_lines[:,0] == float(number_diabats))[0])

    final_position = pseudoH_lines[final_position_index, :2]

    vector_population_array = np.zeros((number_timesteps, number_diabats))
       
    timestep = 0
    initial_index = 0

    for index in range(1, file_length):

        if np.sum(final_position) == np.sum(pseudoH_lines[index][:2]):

            final_index = index + 1

            single_H = pseudoH_lines[initial_index:final_index]

            hamiltonian = np.zeros((number_diabats, number_diabats))

            for element in single_H:

                row = int(element[0]) - 1
                column = int(element[1]) - 1
                quantity = float(element[2])

                if row == column:
                    quantity = quantity/2

                hamiltonian[row,column] = quantity

            lowest_site_energy = np.min(np.diag(hamiltonian))
            hamiltonian = hamiltonian - np.identity(number_diabats)*lowest_site_energy

            hamiltonian = hamiltonian + hamiltonian.T
 
            eigenvals, eigenvecs = np.linalg.eigh(hamiltonian)
            idx = eigenvals.argsort()

            eigenvals = eigenvals[idx]
            eigenvecs = eigenvecs[:, idx]

            active_state = int(active_state_array[timestep]-1)

            active_vector = eigenvecs[:,active_state]
            vector_population_array[timestep, :] = active_vector**2

            timestep += 1
            initial_index = final_index

    return vector_population_array


def get_location_array(population_array, exciton_coms, acceptor_coms, donor_coms, number_diabats, number_excitons):

    number_timesteps = len(population_array)
    location_array = np.zeros((number_timesteps, 3))
    IPR_array = np.zeros((number_timesteps, 3))

    for timestep in range(number_timesteps):

        exciton_location, electron_location, hole_location = diabat_locations(donor_coms, acceptor_coms, acceptor_coms, population_array[timestep], number_diabats, number_excitons)

        exciton_IPR, electron_IPR, hole_IPR = IPR(population_array[timestep], number_diabats, number_excitons)

        location_array[timestep, 0] = exciton_location
        location_array[timestep, 1] = electron_location
        location_array[timestep, 2] = hole_location

        IPR_array[timestep, 0] = exciton_IPR
        IPR_array[timestep, 1] = electron_IPR
        IPR_array[timestep, 2] = hole_IPR

    return location_array, IPR_array

active_state_array = sh_log_parser(f'{target_directory}/run-fssh-{trajectory}/run-sh-1.log')

pseudoH_lines = xyz_parser(f'{target_directory}/run-fssh-{trajectory}/run-pseudo-hamilt-1.xyz', 3)

vector_float = np.vectorize(float)
pseudoH_lines = vector_float(pseudoH_lines)

vector_population_array = inner_loop(pseudoH_lines, active_state_array, len(active_state_array), number_diabats)

location_array, IPR_array = get_location_array(vector_population_array, trajectory_acceptor_xCOMS[int(trajectory)], trajectory_acceptor_xCOMS[int(trajectory)], trajectory_donor_xCOMS[int(trajectory)], number_diabats, number_excitons)


data_file = open(f'{target_directory}/physopt_eigenstate_exciton_locations.txt', 'a')
data_file.write(trajectory + ', ' + str(list(location_array[:, 0]))[1:-1] + '\n')
data_file.close()

data_file = open(f'{target_directory}/physopt_eigenstate_electron_locations.txt', 'a')
data_file.write(trajectory + ', ' + str(list(location_array[:, 1]))[1:-1] + '\n')
data_file.close()

data_file = open(f'{target_directory}/physopt_eigenstate_hole_locations.txt', 'a')
data_file.write(trajectory + ', ' + str(list(location_array[:, 2]))[1:-1] + '\n')
data_file.close()

data_file = open(f'{target_directory}/physopt_eigenstate_exciton_IPR.txt', 'a')
data_file.write(trajectory + ', ' + str(list(IPR_array[:, 0]))[1:-1] + '\n')
data_file.close()

data_file = open(f'{target_directory}/physopt_eigenstate_electron_IPR.txt', 'a')
data_file.write(trajectory + ', ' + str(list(IPR_array[:, 1]))[1:-1] + '\n')
data_file.close()

data_file = open(f'{target_directory}/physopt_eigenstate_hole_IPR.txt', 'a')
data_file.write(trajectory + ', ' + str(list(IPR_array[:, 2]))[1:-1] + '\n')
data_file.close()
