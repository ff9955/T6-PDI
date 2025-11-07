import numpy as np
from file_parsers import *
from xsh_analysis_functions import *
import sys

target_directory = sys.argv[1] 
number_diabats = 12

trajectory_list = [t for t in range(0,1)]
site_energy_array = np.zeros((2000, number_diabats), dtype=float)

for t_count, trajectory in enumerate(trajectory_list):
    print(trajectory)

    diabat_pops = diabatic_populations(f'{target_directory}/run-fssh-{trajectory}/run-coeff-1.xyz', number_diabats)
#    site_energy_lines = xyz_parser(f'{target_directory}/run-fssh-{trajectory}/run-site-1.xyz', 2)

    np.savetxt(f'{target_directory}/3x3_fi_diabat_populations.txt', diabat_pops)

    counter = 0

#    while counter < len(site_energy_lines):
#
#        site_energy_array[counter//number_diabats-1,:] = site_energy_lines[counter : counter + number_diabats ,1]
#        counter += number_diabats
#
#    pseudoH_lines = xyz_parse_first_section(f'{target_directory}/run-fssh-{trajectory}/run-pseudo-hamilt-1.xyz', 3)
#    initial_H = build_sim_H(pseudoH_lines, number_diabats)
#    
#    eigenvals, eigenvecs = get_eigen(initial_H)
#
#    active_state_array = sh_log_parser(f'{target_directory}/run-fssh-{trajectory}/run-sh-1.log')
#    initial_active_state = int(active_state_array[0] - 1)
#
#    print(eigenvecs[:, initial_active_state])
#
#np.savetxt(f'{target_directory}/2phase_site_energies.txt', site_energy_array)
