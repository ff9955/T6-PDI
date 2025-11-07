import numpy as np
from file_parsers import *
from xsh_fast_analysis_functions import *
import sys
import time

target_directory = sys.argv[3] 

trajectory_list = [t for t in range(int(sys.argv[1]),int(sys.argv[2]))]

#bad_trajectories = [46]
#
#for element in bad_trajectories:
#    trajectory_list.remove(element)

CSS_array = np.zeros((len(trajectory_list), 10000))
XT_array = np.zeros((len(trajectory_list), 10000))
INT_array = np.zeros((len(trajectory_list), 10000))

#exciton_IPR_array = np.zeros((len(trajectory_list), 22000))
#hole_IPR_array = np.zeros((len(trajectory_list),  22000))
#electron_IPR_array = np.zeros((len(trajectory_list),  22000))

#exciton_centre_array = np.zeros((len(trajectory_list), 22000))
#hole_centre_array = np.zeros((len(trajectory_list), 22000))
#electron_centre_array = np.zeros((len(trajectory_list),  22000))

#initial_T6_distances = np.loadtxt(f'{target_directory}/T6_initial_xCOMS.txt')
#initial_PDI_distances = np.loadtxt(f'{target_directory}/PDI_initial_xCOMS.txt')

number_diabats = 1300
number_excitons = 50

vector_float = np.vectorize(float)

for t_count, trajectory in enumerate(trajectory_list):
    print(trajectory)

    coulomb_array = xyz_parser(f'{target_directory}/FRZ_DIAGONALS_COULOMB.include', 3)
    coulomb_array = vector_float(coulomb_array[:,-1])

    distance_array = xyz_parser(f'{target_directory}/CT_DISTANCES.include', 3)
    distance_array = vector_float(distance_array[:,-1])
 
    with open(f'{target_directory}/run-fssh-{trajectory}/run-coeff-1.xyz') as coeff_file:

        diabat_counter = 0
        timestep = 0
        coeff_list = []
        population_array = np.zeros((1,number_diabats), dtype=float)

        for line in coeff_file:
            split_line = line.split()

            coeff_list.append(split_line[2:])
            diabat_counter += 1

            if diabat_counter == (number_diabats+2):
                coeff_list = coeff_list[2:]
                coeff_array = np.array(coeff_list, dtype=float)

                coeff_array = np.sum(coeff_array**2, axis=1)
                population_array[0,:] = coeff_array

                CSS, XT, INT = state_type_populations(population_array, distance_array, coulomb_array, 1.0, 1.0, number_diabats)
#                exciton_IPR, electron_IPR, hole_IPR = IPR(population_array[0], number_diabats, number_excitons)
#                exciton_loc, electron_loc, hole_loc = diabat_locations(initial_T6_distances[trajectory], initial_PDI_distances[trajectory], initial_PDI_distances[trajectory], population_array[0], number_diabats, number_excitons)

                CSS_array[t_count,timestep] = CSS
                XT_array[t_count,timestep] = XT
                INT_array[t_count,timestep] = INT

#                exciton_IPR_array[t_count,timestep] = exciton_IPR
#                electron_IPR_array[t_count,timestep] = electron_IPR
#                hole_IPR_array[t_count,timestep] = hole_IPR

#                exciton_centre_array[t_count,timestep] = exciton_loc
#                electron_centre_array[t_count,timestep] = electron_loc
#                hole_centre_array[t_count,timestep] = hole_loc

                timestep += 1
                coeff_list = []
                diabat_counter = 0

zero_row_locs = np.where(XT_array + INT_array == 0)[1]
earliest_zero_index = np.min(zero_row_locs)
print('earliest zero index:', earliest_zero_index)

CSS_array = CSS_array[:, :earliest_zero_index]
XT_array = XT_array[:, :earliest_zero_index]
INT_array = INT_array[:, :earliest_zero_index]

#exciton_IPR_array = exciton_IPR_array[:, :earliest_zero_index]
#electron_IPR_array = electron_IPR_array[:, :earliest_zero_index]
#hole_IPR_array = hole_IPR_array[:, :earliest_zero_index]

#exciton_centre_array = exciton_centre_array[:, :earliest_zero_index]
#electron_centre_array = electron_centre_array[:, :earliest_zero_index]
#hole_centre_array = hole_centre_array[:, :earliest_zero_index]

np.savetxt(f'{target_directory}/0.01fs_CSS_individual_populations.txt', CSS_array)
np.savetxt(f'{target_directory}/0.01fs_XT_individual_populations.txt', XT_array)
np.savetxt(f'{target_directory}/0.01fs_INT_individual_populations.txt', INT_array)

#np.savetxt(f'{target_directory}/e5_2xCT_individual_exciton_IPR.txt', exciton_IPR_array)
#np.savetxt(f'{target_directory}/e5_2xCT_individual_electron_IPR.txt', electron_IPR_array)
#np.savetxt(f'{target_directory}/e5_2xCT_individual_hole_IPR.txt', hole_IPR_array)

#np.savetxt(f'{target_directory}/e5_2xCT_individual_exciton_locations.txt', exciton_centre_array)
#np.savetxt(f'{target_directory}/e5_2xCT_individual_electron_locations.txt', electron_centre_array)
#np.savetxt(f'{target_directory}/e5_2xCT_individual_hole_locations.txt', hole_centre_array)
