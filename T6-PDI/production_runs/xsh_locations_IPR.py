import numpy as np
from file_parsers import *
from xsh_analysis_functions import *
import sys

target_directory = sys.argv[1]

trajectory_list = [t for t in range(0,500)]

#bad_trajectories = [71,197,307,408,475]
#
#for element in bad_trajectories:
#    trajectory_list.remove(element)

exciton_IPR_array = np.zeros((len(trajectory_list), 41000))
donor_IPR_array = np.zeros((len(trajectory_list), 41000))
acceptor_IPR_array = np.zeros((len(trajectory_list), 41000))

exciton_centre_array = np.zeros((len(trajectory_list), 41000))
donor_centre_array = np.zeros((len(trajectory_list), 41000))
acceptor_centre_array = np.zeros((len(trajectory_list), 41000))

trajectory_donor_xCOMS = np.loadtxt(f'{target_directory}/T6_initial_xCOMS.txt')
trajectory_acceptor_xCOMS = np.loadtxt(f'{target_directory}/PDI_initial_xCOMS.txt')

for index, trajectory in enumerate(trajectory_list):
    print(trajectory)

    IPR_array = single_trajectory_IPR(f'{target_directory}/run-fssh-{trajectory}/run-coeff-1.xyz', 60, 10, 5, 10)
    centre_array = single_trajectory_diabat_locations(trajectory_donor_xCOMS[trajectory], trajectory_acceptor_xCOMS[trajectory], trajectory_acceptor_xCOMS[trajectory], f'{target_directory}/run-fssh-{trajectory}/run-coeff-1.xyz', 60, 10, 5, 10)

    IPR_array_length = len(IPR_array[:,0])
    centre_array_length = len(centre_array[:,0])

    exciton_IPR_array[index, :IPR_array_length] = IPR_array[:,0]
    acceptor_IPR_array[index, :IPR_array_length] = IPR_array[:,1]
    donor_IPR_array[index, :IPR_array_length] = IPR_array[:,2]

    exciton_centre_array[index, :centre_array_length] = centre_array[:,0]
    acceptor_centre_array[index, :centre_array_length] = centre_array[:,1]
    donor_centre_array[index, :centre_array_length] = centre_array[:,2]

zero_row_locs = np.where((exciton_IPR_array) == 0)[1]
IPR_zero_index = np.min(zero_row_locs)

zero_row_locs = np.where(exciton_centre_array == 0)[1]
centre_zero_index = np.min(zero_row_locs)

if centre_zero_index < IPR_zero_index:
    earliest_zero_index = centre_zero_index
elif IPR_zero_index < centre_zero_index:
    earliest_zero_index = IPR_zero_index
else:
    earliest_zero_index = IPR_zero_index

exciton_IPR_array = exciton_IPR_array[:, :earliest_zero_index]
donor_IPR_array = donor_IPR_array[:, :earliest_zero_index]
acceptor_IPR_array = acceptor_IPR_array[:, :earliest_zero_index]

exciton_centre_array = exciton_centre_array[:, :earliest_zero_index]
donor_centre_array = donor_centre_array[:, :earliest_zero_index]
acceptor_centre_array = acceptor_centre_array[:, :earliest_zero_index]

np.savetxt(f'{target_directory}/physopt_shortchain_Eb190_individual_exciton_IPR.txt', exciton_IPR_array)
np.savetxt(f'{target_directory}/physopt_shortchain_Eb190_individual_electron_IPR.txt', acceptor_IPR_array)
np.savetxt(f'{target_directory}/physopt_shortchain_Eb190_individual_hole_IPR.txt', donor_IPR_array)

np.savetxt(f'{target_directory}/physopt_shortchain_Eb190_individual_exciton_locations.txt', exciton_centre_array)
np.savetxt(f'{target_directory}/physopt_shortchain_Eb190_individual_electron_locations.txt', acceptor_centre_array)
np.savetxt(f'{target_directory}/physopt_shortchain_Eb190_individual_hole_locations.txt', donor_centre_array)

#TESTING THE IPR FUNCTION BY SEEING IF IT GIVES THE CORRECT IPR, GIVEN SOME SIMPLE POPULATION DISTRIBUTIONS
#----------------------------------------------------------------------------------------------------------
#artificial_frame = np.zeros(420)
#artificial_frame[np.arange(0,420,20)] = 1/21
#artificial_frame[-1] = 1/21
#
#exciton_ipr, acceptor_ipr, donor_ipr = IPR(artificial_frame, 420, 20)
#print(exciton_ipr, acceptor_ipr, donor_ipr)


