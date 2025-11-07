import numpy as np
import sys

simulation = sys.argv[1]
file_names = [f'{simulation}_active_rel_energy-noHband.txt']

bad_trajectories = [19, 32, 45, 103, 226, 340, 399]

for file in file_names:

    data_list = []

    with open(file) as full_file:

        all_lines = full_file.readlines()
        all_lines = list(map(lambda x: x.split(', '), all_lines))

        for line in all_lines:
            if int(line[0]) in bad_trajectories:
                continue
            else:
                data_list.append(line)

    array_lengths = list(map(lambda x: len(x), data_list))
    minimum_length = 20001

    data_list = list(map(lambda x: x[:minimum_length], data_list))
    data_array = np.array(data_list, dtype=float)

    sorted_traj_indices = data_array[:,0].argsort()
    data_array = data_array[sorted_traj_indices, :]

    np.savetxt(file[:-4] + '-p.txt', data_array)
