import numpy as np
from xsh_analysis_functions import get_general_output
from file_parsers import sh_log_parser
import sys

target_directory = sys.argv[1]
trajectory_list = [t for t in range(0,60)]

number_diabats = 12
number_timesteps = 8000

surface_array = np.zeros(number_diabats)
adiabat_bw_array = np.zeros(number_diabats)

temperature = 300
kB = 3.1668114E-6

for t_count, trajectory in enumerate(trajectory_list):
    print(trajectory)

    active_state_array = sh_log_parser(f'{target_directory}/run-fssh-{trajectory}/run-sh-1.log')
    active_state_array = active_state_array[:number_timesteps]

    for timestep in range(len(active_state_array)):

        active_state = int(active_state_array[timestep] - 1)
        surface_array[active_state] = surface_array[active_state] + 1

    adiabat_energy = get_general_output(f'{target_directory}/run-fssh-{trajectory}/run-adiab-1.xyz', 2, 1, 12)
    adiabat_energy = adiabat_energy[:number_timesteps, :]
    
    for k in range(len(adiabat_energy)):

        adiabat_energy[k] = adiabat_energy[k] - adiabat_energy[k,0]

    adiabat_bw = np.exp(-adiabat_energy/(temperature*kB))
    adiabat_bw = np.mean(adiabat_bw, axis=0)
    adiabat_bw_array = adiabat_bw_array + adiabat_bw

adiabat_bw_array = adiabat_bw_array/len(trajectory_list)
adiabat_bw_array = adiabat_bw_array/np.sum(adiabat_bw_array)

surface_array = surface_array/np.sum(surface_array)

#print(adiabat_bw_array)
print(surface_array)
