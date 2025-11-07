import numpy as np
from file_parsers import *
import sys
from numba import jit, int64, float64, types
from numba.experimental import jitclass
from numba.typed import List

target_directory = sys.argv[1]
trajectory_list = [sys.argv[2]]

spec = [('eigenvectors', float64[:,:]), ('eigenvalues', float64[:]), ('active_state', int64), ('number_excitons', int64), ('iCT_index', int64)]

@jitclass(spec)
class surface_wavefunctions():

    def __init__(self, eigenvectors, eigenvalues, active_state, number_excitons, iCT_index):

        self.eigenvectors = eigenvectors
        self.eigenvalues = eigenvalues
        self.number_excitons = number_excitons
        self.iCT_index = iCT_index
        self.active_state = active_state
    
    def boltzmann_weighting(self, temperature, energy_difference):
        #calculating Boltzmann factor of eigenstate of H(R) relative to active state

        kB = 3.1668114E-6
        thermal_energy = kB*temperature
        b = np.exp(-energy_difference/thermal_energy)

        return min([1, b]) #state lower in energy --> factor set to 1

    def distance_resolved_nace(self, nace_array, temperature, hole_coords, electron_coords, CT_threshold):

        active_eigenvector = self.eigenvectors[:,self.active_state]
        active_energy = self.eigenvalues[self.active_state]

        eh_distances = np.zeros(len(self.eigenvectors) - self.number_excitons)
        counter = 0

        #setting e-h distance that corresponds to each CT-diabat, indexing goes: hole + all e-sites, 2nd hole + all e-sites etc...
        for h_coord in hole_coords:
            for e_coord in electron_coords:

                diabat_distance = e_coord - h_coord
                eh_distances[counter] = diabat_distance
                counter += 1

        distance_bins = np.arange(0, np.max(eh_distances) + 20, 20)[1:]
        distance_resolved_naces = [List.empty_list(types.float64) for b in distance_bins]
        distance_resolved_weights = [List.empty_list(types.float64) for b in distance_bins]
        distance_resolved_ediffs = [List.empty_list(types.float64) for b in distance_bins]

        for index in range(len(self.eigenvectors)):

            if index != self.active_state: #iterating over all eigenstates except active state

                nace = nace_array[index] #indexing of nace corresponds to indexing of eigenstates in diagonalised H (checked)
                energy = self.eigenvalues[index]
                eigenvec = self.eigenvectors[:,index]
                
                energy_difference = energy - active_energy
                bw = self.boltzmann_weighting(temperature, energy_difference) #bw = boltzmann weight wrt active state

                ct_diabat_pops = eigenvec[:-self.number_excitons]**2
                ict_population = ct_diabat_pops[360]

                #this is where I specify the required character of the eigenstate whose nace I want to calculate
                if np.sum(ct_diabat_pops) < CT_threshold: #criterion for full CT-states
                    continue
#                if (np.sum(ct_diabat_pops) > CT_threshold) or (np.sum(ct_diabat_pops) < (1-CT_threshold)):
#                    continue
                else:
                    ct_diabat_pops = ct_diabat_pops/np.sum(ct_diabat_pops)

                weighted_distances = ct_diabat_pops*eh_distances
                avg_eh_distance = np.sum(weighted_distances) #expectation value of eigenstate's e-h distance

                #checking which e-h distance bin this eigenstate falls into
                for j in range(len(distance_bins)):
                    distance_bin = distance_bins[j]

                    if avg_eh_distance < distance_bin:
                        distance_resolved_naces[j].append(abs(nace))
                        distance_resolved_weights[j].append(bw)
                        distance_resolved_ediffs[j].append(abs(energy_difference))
                        #append observables to list (within list) corresponding to bin

                        break

        return distance_resolved_naces, distance_resolved_weights, distance_resolved_ediffs

    def IPR_resolved_nace(self, nace_array, temperature, CT_threshold):
        #basically the same function as above, but distance bins replaced with exciton IPR bins

        active_eigenvector = self.eigenvectors[:,self.active_state]
        active_energy = self.eigenvalues[self.active_state]

        IPR_bins = np.arange(0, 12, 2)[1:]
        eIPR_resolved_naces = [List.empty_list(types.float64) for b in IPR_bins]
        eIPR_resolved_weights = [List.empty_list(types.float64) for b in IPR_bins]
        eIPR_resolved_ediffs = [List.empty_list(types.float64) for b in IPR_bins]

        for index in range(len(self.eigenvectors)):

            if index != self.active_state:

                nace = nace_array[index]
                energy = self.eigenvalues[index]
                eigenvec = self.eigenvectors[:,index]
                
                energy_difference = energy - active_energy
                bw = self.boltzmann_weighting(temperature, energy_difference)

                ct_diabat_pops = eigenvec[:-self.number_excitons]**2

                #this is where I specify the required character of the eigenstate whose nace I want to calculate
#                if np.sum(ct_diabat_pops) < CT_threshold: #criterion for full CT-states
#                    continue
                if (np.sum(ct_diabat_pops) > CT_threshold) or (np.sum(ct_diabat_pops) < (1-CT_threshold)):
                    continue
                else:
                    unnormalised_exciton_populations = eigenvec[-self.number_excitons:]**2
                    normalised_exciton_populations = unnormalised_exciton_populations/np.sum(unnormalised_exciton_populations)

                exciton_IPR = 1/np.sum(normalised_exciton_populations**2)

                for j in range(len(IPR_bins)):
                    IPR_bin = IPR_bins[j]

                    if exciton_IPR < IPR_bin:

                        eIPR_resolved_naces[j].append(abs(nace))
                        eIPR_resolved_weights[j].append(bw)
                        eIPR_resolved_ediffs[j].append(abs(energy_difference))
                        break

        return eIPR_resolved_naces, eIPR_resolved_weights, eIPR_resolved_ediffs

    def get_hybrid_info(self, hole_coords, electron_coords, CT_threshold, active_state=False):
        #function for extracting properties from all hybrid eigenstates of H(R), active_state keyword determines whether you want to sample just active hybrid states or all hybrid states

        hybrid_info = List()
        eh_distances = np.zeros(len(self.eigenvectors) - self.number_excitons)

        number_CT_states = len(eh_distances)
        counter = 0
        number_DA_molecules = len(hole_coords)

        #getting e-h distances of CT-states as before
        for h_coord in hole_coords:
            for e_coord in electron_coords:

                diabat_distance = e_coord - h_coord
                eh_distances[counter] = diabat_distance
                counter += 1

        lowest_energy = np.min(self.eigenvalues)
        #ground state energy

        for index in range(len(self.eigenvectors)):

            if (active_state == True) and (index != self.active_state):
                continue

            energy = self.eigenvalues[index]
            eigenvec = self.eigenvectors[:,index]
            
            energy_difference = energy - lowest_energy
            #hybrid state energy wrt ground state

            ct_diabat_pops = eigenvec[:-self.number_excitons]**2

            #eigenstate is hybrid if 0.05 < CT_pop < 0.95
            if (np.sum(ct_diabat_pops) > CT_threshold) or (np.sum(ct_diabat_pops) < (1-CT_threshold)):
                continue
            else:
                unnormalised_exciton_populations = eigenvec[-self.number_excitons:]**2
                normalised_exciton_populations = unnormalised_exciton_populations/np.sum(unnormalised_exciton_populations)
                ct_diabat_pops = ct_diabat_pops/np.sum(ct_diabat_pops)

            exciton_IPR = 1/np.sum(normalised_exciton_populations**2)
            exciton_distance = np.sum(normalised_exciton_populations*electron_coords) - 92
            #exciton distance assumes that exciton is on acceptor phase, need to swap 'electron_coords' if no longer true

            weighted_distances = ct_diabat_pops*eh_distances
            avg_eh_distance = np.sum(weighted_distances)

            #calculating carrier IPR of hybrid eigenstate
            acceptor_site_pops = np.zeros(number_DA_molecules)
            donor_site_pops = np.zeros(number_DA_molecules)

            for j in range(number_DA_molecules):

                donor_sites = eigenvec[j*number_DA_molecules : j*number_DA_molecules + number_DA_molecules]**2
                #get the expansion coeffs of all diabats where the donor's carrier is on the same donor molecule

                acceptor_sites = eigenvec[np.arange(j, number_CT_states ,number_DA_molecules)]**2
                #get the expansion coeffs of all diabats where the acceptor's carrier is on the same molecule

                acceptor_population = np.sum(acceptor_sites)
                donor_population = np.sum(donor_sites)
                #get these selected diabats' populations

                donor_site_pops[j] = donor_population
                acceptor_site_pops[j] = acceptor_population
                #then append the new states' populations to an external array

            acceptor_site_pops = acceptor_site_pops/np.sum(acceptor_site_pops)
            donor_site_pops = donor_site_pops/np.sum(donor_site_pops)
            #re-normalise these populations

            acceptor_IPR = 1/np.sum(acceptor_site_pops**2)
            donor_IPR = 1/np.sum(donor_site_pops**2)
            #and then calculate their IPR, as if they are the only states that make up the wavefunction

            hybrid_info.append(List([exciton_IPR, acceptor_IPR, donor_IPR, exciton_distance, avg_eh_distance, energy_difference]))
            #append list of single state's properties to list of lists

            return hybrid_info
 

@jit
def build_diagonalise(pseudoH_lines, number_diabats):
    #after single Hamiltonian of run-pseudo-hamilt-1.xyz file is parsed, pass the array of strings into this function + the size of the Hilbert space to generate the full H(R) as a 2D numpy array and then diagonalise it

    hamiltonian = np.zeros((number_diabats, number_diabats))

    for element in pseudoH_lines:

       row = int(element[0]) - 1
       column = int(element[1]) - 1
       quantity = element[2]

       if row == column:
           quantity = quantity/2

       hamiltonian[row,column] = quantity

    hamiltonian = hamiltonian - np.identity(number_diabats)*hamiltonian[360,360]

    hamiltonian = hamiltonian + hamiltonian.T

    eigenvals, eigenvecs = np.linalg.eigh(hamiltonian)
    idx = eigenvals.argsort()

    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]

    return eigenvals, eigenvecs


def inner_loop_DI_resolved(pseudoH_lines, nace_array, number_diabats, number_excitons, hole_coords, electron_coords):

    eigenvals, eigenvecs = build_diagonalise(pseudoH_lines, number_diabats)
    #build and diagonalise H(R) parsed from single frame
     
    active_state = int(nace_array[0]-1)
    active_eigenvector = eigenvecs[:,active_state]

    if (np.sum(active_eigenvector[-number_excitons:]**2) < 0.95) and (np.sum(active_eigenvector[:-number_excitons]**2) > 0.05): #and (np.sum(active_eigenvector[-number_excitons+4:]**2) > 0.80): #criterion for non-interfacial xt-dissoc

        surface_info = surface_wavefunctions(eigenvecs.astype(np.float64), eigenvals.astype(np.float64), np.int64(active_state), np.int64(number_excitons), np.int64(360))
#        r_naces, r_weights, r_ediffs = surface_info.distance_resolved_nace(nace_array[2:].astype(np.float64), np.float64(300), hole_coords.astype(np.float64), electron_coords.astype(np.float64), np.float64(0.95))
        hybrid_info = surface_info.get_hybrid_info(hole_coords, electron_coords, np.float64(0.95), True)
        
        return hybrid_info

    elif np.sum(active_eigenvector[-number_excitons:]**2) < 0.95:
        return (True)

    else:
        return (False)

def process_dr_output(output_tuple):
    #function for doing further operations on the nace functions' outputs and sorting them into a single list

    output_lists = [[], [], [], [], [], []] #one inner list for each observable

    for index in range(len(output_tuple[0])): #remember that nace functions output of list of lists where each list has No. elements = No. bins, we're iterating over each bin here

        if output_tuple[0][index]:
            output_lists[0].append(np.mean(output_tuple[0][index])) #average of observable '0' in bin 'index'
            output_lists[1].append(np.mean(output_tuple[1][index]))
            output_lists[2].append(np.sum(np.array(output_tuple[0][index])*np.array(output_tuple[1][index])))
            #sum of boltzmann-weighted nace (bin 0 x bin 1) in bin 'index'
            output_lists[3].append(np.max(output_tuple[0][index]))
            output_lists[4].append(np.sum(output_tuple[1][index]))
            output_lists[5].append(np.mean(output_tuple[2][index]))

        else:
            output_lists[0].append(None) #appending 'None' if list is empty to prevent inhomogeneity
            output_lists[1].append(None)
            output_lists[2].append(None)
            output_lists[3].append(None)
            output_lists[4].append(None)
            output_lists[5].append(None)

    r_nace_output, r_weight_output, r_wnace_output, r_mnace_output, r_sweight_output, r_ediff_output = output_lists[0], output_lists[1], output_lists[2], output_lists[3], output_lists[4], output_lists[5]

    return r_nace_output, r_weight_output, r_wnace_output, r_mnace_output, r_sweight_output, r_ediff_output

#total_nace_list = []
#total_weight_list = []
#total_ediff_list = []
#total_wnace_list = []
#total_mnace_list = []
#total_weightsum_list = []
total_hybrid_info_list = []

donor_xCOMS = np.loadtxt(f'{target_directory}/T6_initial_xCOMS.txt')
acceptor_xCOMS = np.loadtxt(f'{target_directory}/PDI_initial_xCOMS.txt')

for trajectory in trajectory_list:

    with open(f'{target_directory}/run-fssh-{trajectory}/run-pseudo-hamilt-1.xyz') as H_file:
        with open(f'{target_directory}/run-fssh-{trajectory}/run-nace-active-1.xyz') as nace_file:

            H_file.readline()
            H_file.readline() #skipping first 2 lines of pseudo-hamilt file
            parsed_H = []
            interfacial = False
            counter = 0

            for line in H_file:

                split_line = line.split()
                if len(split_line) == 3:
                    parsed_H.append(split_line)

                else:
                    parsed_H = np.array(parsed_H, dtype=float)

                    if counter == 0: #block to skip first H in physopt sim, bc nace at t=0 isn't printed for this simulation
                        counter += 1

                        parsed_H = []
                        H_file.readline()

                        continue

                    nace_file.readline()
                    nace_file.readline()
                    parsed_nace = nace_file.readline().split()
                    parsed_nace = np.array(parsed_nace, dtype=float)

                    output_tuple = inner_loop_DI_resolved(parsed_H, parsed_nace, 420, 20, donor_xCOMS[int(trajectory)], acceptor_xCOMS[int(trajectory)])

                    #the following conditional block reads to the different possible outputs of the inner loop function, and uses them to decide on which batch of naces to save and which to skip over

                    if type(output_tuple) != bool:
                        interfacial = True

#                        mean_naces, mean_weights, mean_wnaces, max_naces, weight_sums, mean_ediffs = process_dr_output(output_tuple)
                        hybrid_info = output_tuple
                        total_hybrid_info_list.append(hybrid_info)

#                        total_nace_list.append(mean_naces)
#                        total_weight_list.append(mean_weights)
#                        total_ediff_list.append(mean_ediffs)
#                        total_mnace_list.append(max_naces)
#                        total_weightsum_list.append(weight_sums)
#                        total_wnace_list.append(mean_wnaces)

#                        for element in hybrid_info:
#                            total_hybrid_info_list.append(list(element))


                        #with else block commented out, we are saving naces from every timestep which satisfies active state requirements, not just those preceding dissociation/hybridisation of active state

#                    else:
#                        if (output_tuple == True) and (interfacial == True):
#
#                            total_nace_list.append(mean_naces)
#                            total_weight_list.append(mean_weights)
#                            total_ediff_list.append(mean_ediffs)
#                            total_wnace_list.append(mean_wnaces)
#                            total_mnace_list.append(max_naces)
#                            total_weightsum_list.append(weight_sums)
#
#                            for element in hybrid_info:
#                                total_hybrid_info_list.append(list(element))
#
#
#                            interfacial = False
#
#                        else:
#                            interfacial = False

                    parsed_H = []
                    H_file.readline()
                    counter += 1

#observables_list = [total_nace_list, total_weight_list, total_wnace_list, total_mnace_list, total_weightsum_list, total_ediff_list]
#text_file_list = ['Ir_XT-H_naces.txt', 'Ir_XT-H_weights.txt', 'Ir_XT-H_wnaces.txt', 'Ir_XT-H_mnaces.txt', 'Ir_XT-H_weight_sums.txt', 'Ir_XT-H_ediffs.txt']
#
##obs_bins = np.arange(0,220,20)[1:]
#obs_bins = np.arange(0,12,2)[1:]
#
#for j in range(len(observables_list)):
#    observable = observables_list[j]
#
#    for ob in range(len(obs_bins)):
#
#        file = open(f'epsilon5_3xCT_{obs_bins[ob]}_{text_file_list[j]}', 'a')
#
#        for timestep in range(len(observable)):
#
#            file.write(str(observable[timestep][ob]) + '\n')
#
#        file.close()

file = open('physopt_retry_active_hybrid_info.txt', 'a')
for j in range(len(total_hybrid_info_list)):

    if total_hybrid_info_list[j] == None:
        continue

    file.write(str(total_hybrid_info_list[j])[1:-1] + '\n' )

file.close()
