import numpy as np
from file_parsers import *
from xsh_analysis_functions import *
import sys

target_directory = sys.argv[1] 

trajectory_list = [t for t in range(0,500)]

#bad_trajectories = [145]
#
#for element in bad_trajectories:
#    trajectory_list.remove(element)

#CALCULATING DIFFERENT DIABAT POPULATIONS FOR INDIVIDUAL TRAJECTORIES, PASSING THEM INTO BIG ARRAYS WITH DIM: NUMBER_TRAJ X NUMBER_PRINTED_TIMESTEPS

#initial_diabat_pops = np.zeros((len(trajectory_list), 1300))
#CSS_array = np.zeros((len(trajectory_list), 41000))
#XT_array = np.zeros((len(trajectory_list), 41000))
#INT_array = np.zeros((len(trajectory_list), 41000))
##FS_array = np.zeros((len(trajectory_list), 41000))
#
#for t_count, trajectory in enumerate(trajectory_list):
#    print(trajectory)
#
#    diabat_pops = diabatic_populations(f'{target_directory}/run-fssh-{trajectory}/run-coeff-1.xyz', 1300)
#    initial_diabat_pops[t_count,:] = diabat_pops[0,:]
#    CSS, FS, XT, INT = state_type_populations(diabat_pops, f'{target_directory}/CT_DISTANCES.include', f'{target_directory}/FRZ_DIAGONALS_COULOMB.include', 2.67, 2.4, 1300)
#
#    CSS_array[t_count, :len(CSS)] = CSS
#    XT_array[t_count, :len(XT)] = XT
#    INT_array[t_count, :len(INT)] = INT
##    FS_array[t_count, :len(FS)] = FS
#
#avg_diabat_pops = np.mean(initial_diabat_pops, axis=0)
##np.savetxt(f'{target_directory}/e5_2xCT_avg_diabats_init.txt', avg_diabat_pops)
#
##PUT WHATEVER YOU WANT AS NUMBER_TIMESTEPS IN THE ARRAY, I FIND THE SMALLEST UNASSIGNED TIMESTEP AND TRUNCATE ALL TRAJECTORIES HERE
#
#zero_row_locs = np.where((XT_array + INT_array + CSS_array) == 0)[1]
#earliest_zero_index = np.min(zero_row_locs)
##earliest_zero_index = 20000
#print('earliest zero index:', earliest_zero_index)
#
#CSS_array = CSS_array[:, :earliest_zero_index]
#XT_array = XT_array[:, :earliest_zero_index]
#INT_array = INT_array[:, :earliest_zero_index]
##FS_array = FS_array[:, :earliest_zero_index]
#
#np.savetxt(f'{target_directory}/e5_2xCT_2D_0.01fs_CSS_individual_populations.txt', CSS_array)
#np.savetxt(f'{target_directory}/e5_2xCT_2D_0.01fs_XT_individual_populations.txt', XT_array)
#np.savetxt(f'{target_directory}/e5_2xCT_2D_0.01fs_INT_individual_populations.txt', INT_array)
#np.savetxt(f'{target_directory}/physical_diabatic_FS_individual_populations.txt', FS_array)

initial_XT_char_array = trajectory_initial_XT(trajectory_list, 'run-sh-1.log', 'run-pseudo-hamilt-1.xyz', target_directory, 61, 10, recomb=True)
np.savetxt(f'{target_directory}/initial_XT_characters-physopt_shortchain_recomb.txt', initial_XT_char_array)

#THEN EXPRESSING THESE ADIABATS' ENERGIES WRT TO THE NEUTRAL GROUND STATE ENERGY AT THAT PARTICULAR GEOMETRY

neutral_energies = np.loadtxt(f'{target_directory}/total_neutral_energies.txt')
non_covalent_energies = np.loadtxt(f'{target_directory}/non_covalent_energies.txt')
selected_geoms = xyz_parser(f'{target_directory}/phys_selected_geoms_adiabats.txt', 2)

vector_int = np.vectorize(int)
selected_geoms = vector_int(selected_geoms[:,0])

CT_energy = 0 #0.086424 for 0.4eV offset in 1D 6T-PDI #0.06878433208233997 for an offset of 0.87 eV

for index, geom_number in enumerate(selected_geoms):

    neutral_energy = neutral_energies[geom_number]
    non_covalent_energy = non_covalent_energies[geom_number]

    adiabat_energy = initial_XT_char_array[0][index]

    shifted_adiabat_energy = adiabat_energy - (neutral_energy - non_covalent_energy) + CT_energy

    initial_XT_char_array[0][index] = shifted_adiabat_energy

np.savetxt(f'{target_directory}/initial_XT_chars_rel-to-neutral-physopt_shortchain_recomb.txt', initial_XT_char_array)

#what I'm doing below with the reference energies, is reading in the eigenvalues of geom 293's Hamiltonian, which were calculated in get_adiabat_pulse.py with the site energies shifted relative to the neutral ground state, so I can compare these eigenvalues with those directly from the unaltered X-SH Hamiltonian, whose eigenvals are then shifted relative to the neutral ground-state. The result is that they should be equal and the difference I am printing below should therefore be zero. The reason I am doing this is because I want to plot the energy (relative to neutral state) distribution of the initial adiabats, and I want to see if I can just get the adiabats' relative energies by just shifting them directly, instead of shifting each individual site energy. Turns out I can.

#reference_energies = np.loadtxt('geom_293_eigenvals.txt')

#pseudoH_lines = xyz_parse_first_section('{target_directory}/run-fssh-0/run-pseudo-hamilt-1.xyz', 3)
#H = build_sim_H(pseudoH_lines, 420)
#eigenvals, eigenvecs = get_eigen(H)

#eigenvals = eigenvals - (neutral_energies[293] - non_covalent_energies[293]) + CT_energy

#print(abs(eigenvals - reference_energies))

#CALCULATING TIME-DEPENDENT ENERGIES OF ADIABATS INSIDE DIFFERENT TYPES OF EIGENSTATE BANDS

#active_states = sh_log_parser('{target_directory}/run-fssh-9/run-sh-1.log')
#np.savetxt('{target_directory}/traj9/physopt-traj9-active-states.txt', np.array(active_states))
#
#H_list = partition_H('{target_directory}/run-fssh-9/run-pseudo-hamilt-1.xyz')
#simulation_times = [time*0.05*10 for time in range(0,len(H_list))]
#
#H_list = H_list[1200:2000:20]
#simulation_times = simulation_times[1200:2000:20]
#
#for index in range(len(H_list)):
#    np.savetxt(f'{target_directory}/traj9/physopt-H-trajectory0-{simulation_times[index]}fs', H_list[index])

#CSS_bands, XT_bands, CT_minus_CSS_bands, hybrid_bands = trajectory_band_structures('run-fssh-399/run-pseudo-hamilt-1.xyz', 'run-fssh-399/CT_DISTANCES.include', 'run-fssh-399/FRZ_DIAGONALS_COULOMB.include', 0.8)
#
#band_list = [CSS_bands, XT_bands, CT_minus_CSS_bands, hybrid_bands]
#band_names = ['phys_CSS_bands_399', 'phys_XT_bands_399', 'phys_CT_minus_CSS_bands_399', 'phys_hybrid_bands_399']

#THEN ONLY TAKING THE VALUES OF THE BAND EDGES AT EACH PRINTED TIMESTEP, PASSING THESE INTO 2-ROW, N_TIMESTEP-COLUMN ARRAYS, THEN SAVING THEM

#for index in range(len(band_list)):
#
#    chosen_band = band_list[index]
#    boundary_array = np.zeros((2,len(chosen_band)))
#
#    for index2 in range(len(chosen_band)):
#        
#        boundary_array[0,index2] = chosen_band[index2][0]
#        boundary_array[1,index2] = chosen_band[index2][-1]
#
#    np.savetxt(f'{band_names[index]}.txt', boundary_array)

#CALCULATING INITIAL IPR VALUES OF TRAJECTORIES, SAVING RESULTING ARRAY OF EXCITON, HOLE AND ELECTRON IPR VALS

#trajectory_donor_xCOMS = np.loadtxt(f'{target_directory}/T6_initial_xCOMS.txt')
#trajectory_acceptor_xCOMS = np.loadtxt(f'{target_directory}/PDI_initial_xCOMS.txt')
#
#initial_IPR_array = trajectory_initial_IPR(trajectory_list, 'run-sh-1.log', 'run-pseudo-hamilt-1.xyz', target_directory, 420, 20)
#np.savetxt(f'{target_directory}/e3.5_offset0.87_initial_IPR.txt', initial_IPR_array)
#
#charge_carrier_centres = np.zeros((len(trajectory_list), 3))
#    
#for index, trajectory in enumerate(trajectory_list):
#
#    exciton_centre, acceptor_centre, donor_centre = initial_diabat_locations(trajectory_donor_xCOMS[trajectory], trajectory_acceptor_xCOMS[trajectory], trajectory_acceptor_xCOMS[trajectory], f'{target_directory}/run-fssh-{trajectory}/run-sh-1.log', f'{target_directory}/run-fssh-{trajectory}/run-pseudo-hamilt-1.xyz', 420, 20)
#
#    charge_carrier_centres[index, 0] = exciton_centre
#    charge_carrier_centres[index, 1] = acceptor_centre
#    charge_carrier_centres[index, 2] = donor_centre
#
#np.savetxt(f'{target_directory}/e3.5_offset0.87_carrier_centres.txt', charge_carrier_centres)
