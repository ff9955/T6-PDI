import numpy as np
from xsh_analysis_functions import get_general_output
from file_parsers import sh_log_parser
import sys

target_directory = sys.argv[1]
trajectory_list = [t for t in range(199,200)]

#bad_trajectories = [19, 32, 45]
#
#for element in bad_trajectories:
#    trajectory_list.remove(element)

number_diabats = 420
number_excitons = 0
number_timesteps = 19500

surface_array = np.zeros((number_timesteps, number_diabats))
adiabat_pop_array = np.zeros((number_timesteps, number_diabats))

for t_count, trajectory in enumerate(trajectory_list):
    print(trajectory)

    active_state_array = sh_log_parser(f'{target_directory}/run-fssh-{trajectory}/run-sh-1.log')
    active_state_array = active_state_array[:number_timesteps]

    for timestep in range(len(active_state_array)):

        active_state = int(active_state_array[timestep] - 1)
        surface_array[timestep, active_state] = surface_array[timestep, active_state] + 1

    adiabatic_pops = get_general_output(f'{target_directory}/run-fssh-{trajectory}/run-adiab_pop-1.xyz', 3, 2, 420)
    adiabatic_pops = adiabatic_pops[:number_timesteps, :]

    adiabat_pop_array[:,:] = adiabat_pop_array[:,:] + adiabatic_pops

for timestep in range(len(surface_array)):

    surface_array[timestep] = surface_array[timestep]/np.sum(surface_array[timestep])

adiabat_pop_array = adiabat_pop_array/len(trajectory_list)

squared_deviation = (surface_array - adiabat_pop_array)**2
mean_squared_deviation = np.mean(squared_deviation, axis=0)

RMSE = np.sqrt(mean_squared_deviation)

mean_surface_populations = np.mean(surface_array, axis=0)
mean_adiabat_populations = np.mean(adiabat_pop_array, axis=0)

np.savetxt(f'{target_directory}/e5_2xCT_mean_adiabat_populations.txt', mean_adiabat_populations)
np.savetxt(f'{target_directory}/e5_2xCT_mean_surface_populations.txt', mean_surface_populations)

#first10_RMSE = RMSE[:10]/mean_surface_populations[:10]
#np.savetxt('physopt_wf_PES_RMSE.txt', first10_RMSE)
