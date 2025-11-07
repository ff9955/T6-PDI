#python script containing a few function you can use to calculate the energy drift of a selected number of trajectories

import numpy as np

def xyz_parser(file_path, number_elements):
    '''
    reads .xyz or .txt file into an array of strings to be processed by another function, the line_elements argument refers to the number of elements your desired lines have, in case different lines in the file have different numebrs of elements
    '''

    with open(file_path) as full_file:
        all_lines = full_file.readlines()

        split_lines = list(map(lambda x: x.split(), all_lines))
        chosen_lines = list(filter(lambda x: True if len(x) == number_elements else False, split_lines))

        filtered_array = np.array(chosen_lines)

    return filtered_array

def energy_drift_gradient(active_directory, trajectory_range, active_atoms):
    '''
    function for calculating the average energy drift of a selection of trajectories in a directory, by calculating the mean of the absolute gradients in energy between all timesteps in each trajectory

    Input: active_directory = string, path to directory containing trajectories; trajectory_range = 2-element iterable with upper and lower bound of trajectories you would like to iterate over; int, active_atoms = number of atoms in the electronically active region

    Output: float, average energy drift
    '''

    trajectory_list = [t for t in range(trajectory_range[0], trajectory_range[1])]
    all_energy_diffs = []

    for trajectory in trajectory_list:

        try:
            ener_lines = xyz_parser(f'{active_directory}/run-fssh-{trajectory}/run-1.ener', 7)
            vector_float = np.vectorize(float)
            ener_array = vector_float(ener_lines)

            total_energies = ener_array[:, 5]

            all_energy_diffs.append(np.diff(total_energies))
            #taking the total energies of all timesteps and calculating the corresponding array of differences in energy between consecutive timesteps

        except:
            print(f'Trajectory {trajectory} has no run-1.ener file')

    timestep_ps = (ener_array[2,1] - ener_array[1,1])/1000
    #getting time in ps between two consecutive timesteps in simulation

    average_energy_diff = 0
    for array in all_energy_diffs:
        average_energy_diff = average_energy_diff + np.mean(abs(array))
        #then average over all energy differences of all trajectories

    average_energy_diff = average_energy_diff/active_atoms
    average_energy_diff = average_energy_diff/timestep_ps
    #divide by number of active atoms and timestep to get the right units

    energy_drift = average_energy_diff/len(all_energy_diffs)
    #doing the energy conservation this way will likely overestimate the energy drift for smaller systems, as you are directly accounting for short-timescale energy fluctuations, which will increase as system size decreases. If you take just the initial and final energies of all the trajectories, these fluctuations will partly cancel eachother out, and reduce the overall energy drift.

    return energy_drift

def energy_drift_diff(active_directory, trajectory_range, active_atoms):
    '''
    function for calculating the energy drift (per atom per ps) by taking the difference between initial and final energy of each trajectory, and averaging over them

    Inputs and outputs are the same as for the 'gradient' function above, the only real difference is the way the energy differences are calculated
    '''

    trajectory_list = [t for t in range(trajectory_range[0], trajectory_range[1])]
    all_energy_diffs = []

    for trajectory in trajectory_list:

        try:
            ener_lines = xyz_parser(f'{active_directory}/run-fssh-{trajectory}/run-1.ener', 7)
            vector_float = np.vectorize(float)
            ener_array = vector_float(ener_lines)
            #parse output from ener file and convert all matrix elements to floats

            total_energies = ener_array[:, 5]
            simulation_time = ener_array[:, 1]
            #index the parsed array such that we now have 1D vectors of time and total energy

            initial_final_E_diff = abs(total_energies[-1] - total_energies[0])
            total_sim_time_ps = abs(simulation_time[-1] - simulation_time[0])/1000
            #calculate difference between intial and final energies and times; convert femtoseconds into picoseconds

            E_diff_per_ps = initial_final_E_diff/total_sim_time_ps
            #then get energy difference (Hartree) per ps

            all_energy_diffs.append(E_diff_per_ps)
            #append to list so we can average over all trajectories' energy differences

        except:
            print(f'Trajectory {trajectory} has no run-1.ener file')

    avg_energy_drift = 0
    for diff in all_energy_diffs:
        avg_energy_drift = avg_energy_drift + diff

    avg_energy_drift = avg_energy_drift/len(all_energy_diffs)
    avg_energy_drift = avg_energy_drift/active_atoms
    #dividng average energy difference by number of atoms in electronically active region to get energy diff per ps per active atom

    return avg_energy_drift

def plot_trajectory_energy(active_directory, trajectory_range):
    '''
    Function for plotting the average total energy of all of the trajectories specified within the chosen range.

    Input: active_directory: string, path to the directory where the trajectories are; trajectory range: two-element iterable, containing the indices of the initial and final trajectories you want to average over

    Output: 2D array, where 1st row contains average energies, and 2nd row contains corresponding smimulation times
    '''

    trajectory_list = [t for t in range(trajectory_range[0], trajectory_range[1])]

    trajectory_energies = []
    trajectory_timesteps = []

    for trajectory in trajectory_list:

        ener_lines = xyz_parser(f'{active_directory}/run-fssh-{trajectory}/run-1.ener', 7)
        vector_float = np.vectorize(float)
        ener_array = vector_float(ener_lines)

        total_energies = ener_array[:, 5]
        simulation_time = ener_array[:, 1]
        #parsing and indexing ener file array to get 1D arrays of time and energy for a single trajectory in the loop

        trajectory_energies.append(total_energies)
        trajectory_timesteps.append(simulation_time)
        #appending these arrays to external lists

    minimum_traj_length = min([len(arr) for arr in trajectory_energies])
    #obtaining the minimum number of timesteps out of all trajectories that have been iterated over

    new_trajectory_energies = np.zeros((len(trajectory_list), minimum_traj_length))
    #defining new energy array where number of columns is the minimum number of timesteps, number of rows is the number of trajectories

    time_energy_array = np.zeros((2, minimum_traj_length))
    time_energy_array[1,:] = trajectory_timesteps[0][:minimum_traj_length]
    #initialising another new array for just time and average energy, assigning top row to time

    for index in range(len(trajectory_energies)):

        new_trajectory_energies[index, :] = trajectory_energies[index][:minimum_traj_length]
        #assigning each trajectory's truncated energy array to the new big energy array

    avg_traj_energies = np.mean(new_trajectory_energies, axis=0)
    time_energy_array[0,:] = avg_traj_energies
    #then calculating the average energy for each timestep, then assiging the resulting average energy array to the new average energy vs time array

    return time_energy_array

def check_timesteps(active_directory, trajectory):
    '''
    Short function that just checks that, for a specific trajectory, no timesteps have been missed out in the run-1.ener file, as this is idicative of memory leaks.
    '''

    try:
        energy_array = xyz_parser(f'{active_directory}/run-fssh-{trajectory}/run-1.ener', 7)
        time_array = energy_array[:,0]

        vector_int = np.vectorize(int)
        time_array = vector_int(time_array)

        for index in range(1,len(time_array)):

            if (time_array[index] - time_array[index-1]) > 1:
                print(time_array[index], trajectory)
                break

    except:
        print(f'encountered binary: traj {trajectory}')


avg_energy_drift = energy_drift_diff('../e5_2xCT_2D_0.03fs', (0,500), 3100) 
print(avg_energy_drift)

trajectory_list = [t for t in range(0,500)]

for trajectory in trajectory_list:
    check_timesteps('../e5_2xCT_2D_0.03fs', trajectory)

