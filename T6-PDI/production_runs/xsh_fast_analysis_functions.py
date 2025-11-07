import numpy as np
from file_parsers import *
from numba import jit

#PYTHON SCRIPT CONTAINING VARIOUS FUNCTIONS FOR ANALYSING X-SH TRAJECTORIES, READ DOCSTRINGS AND COMMENTS TO FIND OUT WHAT EACH FUNCTION DOES AND HOW THEY LINK TOGETHER

def build_sim_H(input_H, dimension):
    '''
    Function that takes the string array of a pseudo-hamiltonian, read-in from the xyz file, and constructs the full Hamiltonian in a numpy array

    Input: input_H: 2D array of string values of one section of the run-pseudo-hamilt-1.xyz file; dimension = int, number of basis states that make up your Hamiltonian

    Output: 2D array of the diabatic X-SH Hamiltonian, where diagonal elements represent site energies and off-diagonals represent electronic couplings between diabats
    '''

    H = np.zeros((dimension, dimension))

    for element in input_H:

        row = int(element[0]) - 1
        #first columns in pseudo-hamilt file denote the index of the row

        column = int(element[1]) - 1
        #second columns denote the index of the column

        quantity = float(element[2])
        #third column therefore contains the value of the matrix element

        if row == column:
            quantity = quantity/2
            #dividing all site energies by 2

        H[row,column] = quantity
        #assiging an element to its corresponding quantity

    H = H + H.T
    #I have only assigned the upper triangle, so I am assigning the rest by adding the transposes together; this works because the Hamiltonian is symmetric, and I halved the site energies earlier

    return H

def get_eigen(matrix):
    '''
    function computes energies and eigenvectors of eigenstates in increasing order of energy

    input: square Hamiltonian matrix

    output: row vector (length D) of eigen energies, matrix (DxD) of orthonormal eigenvectors
    where each column is a single eigenvector
    '''
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    idx = eigenvalues.argsort()

    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    return eigenvalues, eigenvectors

def initial_XT_character(sh_file_path, pseudoH_file_path, number_diabats, number_excitons):
    '''
    Function that calculate the total electronic population of excitonic diabats making up the eigenvectors of the initial surface wavefunction of an X-SH trajectory

    Input: sh_file_path: string, path to run-sh-1.log file; pseudoH_file_path: string, path to run-pseudo-hamilt-1.xyz file; number_diabats: int, number of diabats composing the electronic wavefunction; number_excitons: int, out of the total number of diabats, how many of them are frenkel excitons?

    Output: XT_character: float, sum of populations of the exciton diabats that make up the surface wavefunction
    '''

    sh_lines = xyz_parse_section(sh_file_path, 5, 1, 2)
    initial_adiabat = int(sh_lines[0][-1]) - 1
    #reading first line of the run-sh-1.log file, where it tells you your initial active state

    pseudoH_lines = xyz_parse_first_section(pseudoH_file_path, 3)
    H = build_sim_H(pseudoH_lines, number_diabats)
    #build the initial Hamiltonian

    eigenvals, eigenvecs = get_eigen(H)
    #then diagonalise to get all possible surface wavefunctions

    adiabat_wavefunction = eigenvecs[:,initial_adiabat]
    adiabat_energy = eigenvals[initial_adiabat]
    #the eigenvecs and active state labels are ordered with respect to energy, so you can use the active state label to index the correct eigenvector (i.e surface wavefunction)

    XT_character = np.sum(adiabat_wavefunction[-number_excitons:]**2)
    #calculate total populations of exciton diabats of the active surface wavefunction only

    return XT_character, adiabat_energy

def trajectory_initial_XT(trajectory_list , sh_file_path, pseudoH_file_path, target_dir, number_diabats, number_excitons):
    '''
    Function which calculates the XT-character (shown above) of the initial surface wavefunction of the trajectories specified in the input.
    '''

    XT_character_list = []
    adiabat_energy_list = []

    for trajectory in trajectory_list:

        xt_char, energy = initial_XT_character(f'{target_dir}/run-fssh-{trajectory}/{sh_file_path}', f'{target_dir}/run-fssh-{trajectory}/{pseudoH_file_path}', number_diabats, number_excitons)

        XT_character_list.append(xt_char)
        adiabat_energy_list.append(energy)

    return np.array([adiabat_energy_list, XT_character_list])

@jit
def diabatic_populations(coeff_array, number_diabats):
    '''
    Function that calculates the populations of all diabats that make up the electronic wavefunction throughout a trajectory

    Input: file_path: string, path to coefficient file; number diabats: int, number of diabatic states making up the electronic wavefunction

    Output: 2D array where each row refers to the diabatic populations at a given timestep, and each column refers to the contributing population of a single diabat
    '''

#    coeff_array = xyz_parser(file_path, 4)
#    coeff_array = coeff_array[:, 2:]
#    #reading in all relevant lines of the coefficient file and neglecting everything but the coeff values
#
#    vector_float = np.vectorize(float)
#    coeff_array = vector_float(coeff_array)
#    #convert coeff values to floats

    number_timesteps = len(coeff_array)//number_diabats
    population_array = np.zeros((number_timesteps, number_diabats))
    #initialise array where nrows = number of timesteps, and each row corresponds to the populations of all diabats at a specific timestep

    array_indices = np.arange(0, len(coeff_array), number_diabats)
    #setting the indices in the coeff file so that I'm slicing thorugh diabats of one timestep at a time

    for counter, index in enumerate(array_indices):

        lower_index = index
        upper_index = index + number_diabats
        #iterating over indices used to slice coeff array

        current_coeffs = coeff_array[lower_index:upper_index]
        #getting coeff array section corresponding to one time step

        diabatic_populations = current_coeffs[:,0]**2 + current_coeffs[:,1]**2
        #population is the sum of squared real and imaginary components of the diabat's expansion coefficient

        population_array[counter] = diabatic_populations
        #assign a row of populations to a certain timestep

    return population_array

def adiabatic_populations(file_path, number_diabats):
    '''
    Function that reads the adiabatic populations of the electronic wavefunction, printed out in the run-adiab_pop-1.xyz file, into an array that contains all adiabats' populations for every printed timestep

    Input: file_path: string, path to file that needs to be read, number_diabats: int, number of diabats/adiabats in system

    Output: 2D numpy array where No. columns is number of adiabats, number of rows is number of timesteps, so each row contains the relative populations of all adiabats in the electronic wavefunction for that particular timestep
    '''

    adiabat_array = xyz_parser(file_path, 3)
    adiabat_array = adiabat_array[:, 2]
    #reading file, slicing out population column

    vector_float = np.vectorize(float)
    adiabat_array = vector_float(adiabat_array)
    #converting to floats

    number_timesteps = len(adiabat_array)//number_diabats

    population_array = np.zeros((number_timesteps, number_diabats))

    array_indices = np.arange(0, len(adiabat_array), number_diabats)
    #setting the indices in the adiab file so that I'm slicing through adiabats of one timestep at a time

    for counter, index in enumerate(array_indices):

        lower_index = index
        upper_index = index + number_diabats
        #iterating over indices used to slice coeff array

        current_pops = adiabat_array[lower_index:upper_index]
        #getting adiab array section corresponding to one time step

        population_array[counter] = current_pops
        #assign a row of populations to a certain timestep

    return population_array

@jit
def state_type_populations(diabat_population_array, distance_array, coulomb_array, hIPR, eIPR, number_diabats):
    '''
    Function that divides the diabatic populations of each timestep into distinct types, which can then be individually plotted

    Input: coeff_file_path: string, path to coefficient file; distance_file_path: string, path to trajectory's CT_DISTANCES.include file, showing distances between carriers for different CT-states; coulomb_barrier_file_path: string, path to include file showing coulomb barrier value of each CT-state of trajectory; number_diabats: int, number of diabats making up electronic wavefunction

    Output: three 1D arrays, where each array refers to the combined population of a certain type of diabat as it varies with time
    '''

#    vector_float = np.vectorize(float)
#    coulomb_array = xyz_parser(coulomb_barrier_file_path, 3)
#    coulomb_array = vector_float(coulomb_array[:,-1])
#
#    distance_array = xyz_parser(distance_file_path, 3)
#    distance_array = vector_float(distance_array[:,-1])
    #read the distances and coulomb vals of each CT-state into row vectors

    number_excitons = len(diabat_population_array[0]) - len(coulomb_array)

    exciton_populations = diabat_population_array[:,-number_excitons:]

    exciton_populations = np.sum(exciton_populations, axis=1)
    #calculate the sum of all exciton populations, using the fact that exciton diabats are always the final block in the electronic Hamiltonian

    sorted_coulomb_array = np.sort(coulomb_array)
    interfacial_CT_index = np.where(coulomb_array == sorted_coulomb_array[0])[0]
    #get the index of the interfacial CT-state by getting the index of the coulomb barrier minimum

    interfacial_CT_populations = diabat_population_array[:,interfacial_CT_index]
    #then use this index to get a row vector of interfacial diabat population vs timestep
    interfacial_CT_populations = np.sum(interfacial_CT_populations, axis=1)

    max_coulomb_index = np.where(coulomb_array == max(coulomb_array))[0][0]
    max_coulomb_distance = distance_array[max_coulomb_index]
    #obtain the values of maximum coulomb value and its associated distance

    #thermal_energy =  3.1668114E-6 * 300  #Temperature is HARD CODED TO 300K, YES SIR YES IT IS, CODED BY THE TRENCH BROTHERS

    lattice_distance = distance_array[1] - distance_array[0]
    electron_cloud = (lattice_distance-1)*eIPR #number of lattice distances the particle will spread over when it has average IPR
    hole_cloud = (lattice_distance-1)*hIPR

    #cutoff_distance = np.max(distance_array) - hole_cloud - electron_cloud 
    #defining new minimum eh-distance where e/h have to be at least 2 cloud-distances away from the distance corresponding to the peak of the coulomb potential
    cutoff_distance = max_coulomb_distance

    CSS_indices = np.where(distance_array > cutoff_distance)[0]
   #getting indices of CT-diabats corresponding to CSS, by getting indices of distance larger than new cutoff distance

    #get indices of CT-states that are further apart than that with the maximum coulomb value
    
    #CSS_indices_by_distance = np.where(distance_array > max_coulomb_distance)
    #CSS_indices_by_coulomb = np.where(coulomb_array < (max_coulomb_value - 2*thermal_energy))[0]
    #then get indices of CT-states that have lower coulomb vals (by 2kT) than the maximum coulomb val

    #CSS_indices = np.intersect1d(CSS_indices_by_distance, CSS_indices_by_coulomb)
    #taking the common indices here will give you the indices of CT-states that are to the right of the coulombic maximum, and more stable by 2kT

    CSS_populations = diabat_population_array[:, CSS_indices]
    #get the time-dependent populations of these specific diabats, which we refer to as charge-separated states

    CSS_populations = np.sum(CSS_populations, axis=1)
    #sum them together to get their total diabatic population vs time

    return CSS_populations, exciton_populations, interfacial_CT_populations

def avg_state_populations(trajectory_list, coeff_file_path, distance_file_path, coulomb_barrier_file_path, number_diabats, number_steps):
    '''
    Function that calculates the average diabatic populations which are sepearated by the three types shown in the function above

    Input: trajectory_list: list of trajectories for which you want to calculate the average electronic populations; other inputs are the same as the function above

    Output: three 1D arrays, where the arrays now contain time-dependent diabatic populations, that are averages of the diabatic populations of all trajectories specified in the input
    '''
 
    number_trajectories = len(trajectory_list)

    CSS_array = np.zeros((number_steps, number_trajectories))
    XT_array = np.zeros((number_steps, number_trajectories))
    INT_array = np.zeros((number_steps, number_trajectories))
    A_INT_array = np.zeros((number_steps, number_trajectories))
        #create arrays to hold all trajectories' time-dependent electronic populations, setting all column lengths to the maximum possible number of printed timesteps, which is passed in as an argument


    for index, trajectory in enumerate(trajectory_list):
        print(trajectory)
        #iterating over trajectories specified in the inout list

        diabat_population_array = diabatic_populations(f'run-fssh-{trajectory}/{coeff_file_path}', number_diabats)
        
        CSS, XT, INT, A_INT = state_type_populations(diabat_population_array, f'run-fssh-{trajectory}/{distance_file_path}', f'run-fssh-{trajectory}/{coulomb_barrier_file_path}', number_diabats)
        #calculating electronic populations for a single trajectory

        CSS_array[:len(CSS), index] = CSS
        XT_array[:len(XT), index] = XT
        INT_array[:len(INT), index] = INT
        A_INT_array[:len(A_INT), index] = A_INT
        #after the three types of population vs time have been calculated for this trajectory, insert them into the arrays, where sum trajectories' lengths may be shorter than the maximum possible number of printed timesteps

    #specify an array of the sums of XT and INT populations, so that electronic populations of these arrays are non-zero at all timesteps that have been read in from the simulation, and the only zero array elements are those that weren't assigned

    zero_column_locs = np.where(XT_array == 0)[0]

    if zero_column_locs.size > 0:

        earliest_zero_index = zero_column_locs[0]
        print(earliest_zero_index)
    #find the earliest occurence of an unassigned element in the population array, so that we can remove all trailing zeros, using the shortest-length trajectory to determine where we start to truncate the time-dependent populations

        CSS_array = CSS_array[:earliest_zero_index,:]
        XT_array = XT_array[:earliest_zero_index, :]
        INT_array = INT_array[:earliest_zero_index, :]
        A_INT_array = A_INT_array[:earliest_zero_index, :]

    CSS_average = np.mean(CSS_array, axis=1)
    XT_average = np.mean(XT_array, axis=1)
    INT_average = np.mean(INT_array, axis=1)
    A_INT_average = np.mean(A_INT_array, axis=1)
    #after all arrays have been fully assigned, average across columns to get the average electronic populations

    return CSS_average, XT_average, INT_average, A_INT_average

def partition_H(pseudoH_path):
    '''
    Function that builds the electronic Hamiltonian for every timestep (for which it was printed) in a given trajectory, and then puts them all in a list

    Input: pseudoH_path: string, path to pseudo-hamilt file you want to read in

    Output: pseudoH_list: list of 2D Hamiltonian arrays, with one Hamiltonian for each printed timestep in the trajectory
    '''

    pseudoH_lines = xyz_parser(pseudoH_path, 3)
    #reading in all relevant lines of the pseudo-hamilt file
    vector_float = np.vectorize(float)

    total_pseudoH_array = vector_float(pseudoH_lines)
    initial_position = total_pseudoH_array[0,:2]
    #taking the two indices of the first element of the Hamiltonian, which will correspond to the first diagonal element

    initial_index = 0
    #index from which to begin reading in a section of total_pseudoH_array

    pseudoH_list = []
    for index in range(1, len(total_pseudoH_array)):
        
        if np.sum(initial_position) == np.sum(total_pseudoH_array[index][:2]):
            #conditional statement basically checks if the current line you're on is the first element of a Hamiltonian

            final_index = index
            #if so, specify this as the final index
            
            single_H = total_pseudoH_array[initial_index:final_index]
            #and then read in the hamiltonian that came before the one you've just entered

            pseudoH_list.append(single_H)
            #then append this array slice to the list, as it corresponds to the array of a hamtiltonian at a given timestep

            initial_index = final_index
            #to prepare to read in the next hamiltonian, set the intitial index you'll read in to the final index of the prevous hamiltonian

    return pseudoH_list

def H_band_structure(Hamiltonian, distance_file_path ,coulomb_barrier_file_path ,population_cutoff):
    '''
    Function that diagonalises a Hamiltonian to get the collection of eigenstates, and bins these eigenstates according to the type of diabats that dominate them (or not)

    Input: Hamiltonian: 2D numpy array of electronic Hamiltonian; distance_file_path: string, path to CT_DISTANCES.include file of trajectory you want to read; coulomb_barrier_file_path: string, path to coulomb barrier file of the trajectory; population_cutoff: float, fraction by which a certain type of diabat has to dominate an eigenstate, for you to classify it a certain way.

    Output: 4 lists, containing the energies of the eigenstates that fall into the category represented by that list.
    '''

    eigenvals, eigenvecs = get_eigen(Hamiltonian)
    #diagonalise the input Hamiltonian, to get its eigenvalues and eigenvectors

    eigenvec_transpose = eigenvecs.T

    vector_float = np.vectorize(float)
    coulomb_array = xyz_parser(coulomb_barrier_file_path, 3)
    coulomb_array = vector_float(coulomb_array[:,-1])

    distance_array = xyz_parser(distance_file_path, 3)
    distance_array = vector_float(distance_array[:,-1])
    #read in coulomb values and e-h distances of CT-state diabats into row vectors

    number_excitons = len(eigenvec_transpose[0]) - len(coulomb_array)
    #calculate the number of XT diabats in your system

    XT_bin = []
    CSS_bin = []
    CT_minus_CSS_bin = []
    hybrid_bin = []
    #I assign 4 different types of eigenstate: one dominated by excitons, one by charge-separated CT-states, one by non-separated CT-states, and a hybrid bin which contains states that aren't dominated by any of the diabat types specified above

    for index, vector in enumerate(eigenvec_transpose):

        exciton_populations = vector[-number_excitons:]**2
        exciton_populations = np.sum(exciton_populations)
        #calculate total diabatic population of all excitons

        max_coulomb_index = np.where(coulomb_array == max(coulomb_array))[0][0]
        max_coulomb_value = max(coulomb_array)
        max_coulomb_distance = distance_array[max_coulomb_index]

        thermal_energy =  3.1668114E-6 * 300

        CSS_indices_by_distance = np.where(distance_array > max_coulomb_distance)[0]
        CSS_indices_by_coulomb = np.where(coulomb_array < (max_coulomb_value - 2*thermal_energy))[0]

        CSS_indices = np.intersect1d(CSS_indices_by_distance, CSS_indices_by_coulomb)
        #as before, get the indices of CSS by taking those to the right of the coulomb barrier maximum, and 2kT below it in energy

        CSS_populations = vector[CSS_indices]**2
        CSS_populations = np.sum(CSS_populations)
        #take the CSS diabats and sum their squares to get their total contribution to the electronic population

        CT_populations = vector[:-number_excitons]**2
        CT_populations = np.sum(CT_populations)
        #get the total diabat population of all CT-states - which is every diabat except excitons

        CT_minus_CSS_populations = CT_populations - CSS_populations
        #since all diabat populations should sum to 1, you can get the non-CSS population by subtracting total CT population from the CSS population, which has already been calculated

        if CSS_populations >= population_cutoff:
            CSS_bin.append(eigenvals[index])
        elif exciton_populations >= population_cutoff:
            XT_bin.append(eigenvals[index])
        elif CT_minus_CSS_populations >= population_cutoff:
            CT_minus_CSS_bin.append(eigenvals[index])
        else:
            hybrid_bin.append(eigenvals[index])

        #these conditional statements check to see if the diabat populations of a certain type make up more than a given fraction of the total electronic population, and if they do, they are appended to the corresponding bin

    return CSS_bin, XT_bin, CT_minus_CSS_bin, hybrid_bin

def trajectory_band_structures(pseudoH_path, distance_file_path ,coulomb_barrier_file_path ,population_cutoff):
    '''
    Function that calculates the band structure of the Hamiltonians printed at all timesteps in a given trajectory.

    Input: pseudoH_path: str, path to pseudo-hamilt file, distance_file_path and coulomb_barrier_file_path = str and self-explanatory, population_cutoff: float, fraction of electronic population a type of diabat has to take up for the underlying eigenstate to be assigned as that diabat type.

    Output: 4 different types of lists-of-lists, corresponding to the different diabat types laid out in the previous function, within a list-of-lists, each list contains the energies of eigenstates that correspond to a certain type of diabat at a given timestep
    '''

    H_list = partition_H(pseudoH_path)
    #read all printed Hamiltonians into a list of 2D arrays

    number_diabats = int(H_list[0][-1][0])
    #the total number of diabats in the system will equal the dimension of the Hamiltonian

    CSS_bands = []
    XT_bands = []
    CT_minus_CSS_bands = []
    hybrid_bands = []

    for timestep, H in enumerate(H_list):
        #iterate over all Hamiltonians that were printed in the trajectory

        full_H = build_sim_H(H, number_diabats)
        #build the full Hamiltonian from the 2D array of strings

        CSS_bin, XT_bin, CT_minus_CSS_bin, hybrid_bin = H_band_structure(full_H, distance_file_path, coulomb_barrier_file_path, population_cutoff)
        #bin the Hamiltonian's eigenstates into the different lists

        if CSS_bin: CSS_bands.append(CSS_bin)
        if XT_bin: XT_bands.append(XT_bin)
        if CT_minus_CSS_bin: CT_minus_CSS_bands.append(CT_minus_CSS_bin)
        if hybrid_bin: hybrid_bands.append(hybrid_bin)
        #then append these lists, which are for one timestep only, to a list of band energies for all timesteps

    return CSS_bands, XT_bands, CT_minus_CSS_bands, hybrid_bands

def initial_diabat_projection(pseudoH_file_path, number_diabats, diabat_index):

    '''
    Function that calculates the projection (a.k.a diabatic populations) of a chosen diabat onto all surface wavefunctions of the t=0 Hamiltonian of a single trajectory

    Input: pseudoH_file_path: string, path to pseudo-hamiltonian file; number_diabats: int, total number of diabats in active region; diabat_index: index of diabat whose projection you want (index starts at zero)

    Output: 2D array, top row is adiabat energies, bottom is diabat projection onto these adiabats
    '''

    diabat_projections = np.zeros((2,number_diabats))

    pseudoH_lines = xyz_parse_first_section(pseudoH_file_path, 3)
    H = build_sim_H(pseudoH_lines, number_diabats)
    #build the initial Hamiltonian

    eigenvals, eigenvecs = get_eigen(H)                                                                                                  #then diagonalise to get all possible surface wavefunctions
    
    for index in range(len(eigenvals)):

        vector = eigenvecs[:, index]
        energy = eigenvals[index]
        #iterating through all eigenstates of H, getting eigenvector and eigen energy

        diabat_projection = vector[diabat_index]**2
        #getting diabat projection (population) of the indexed eigenvector

        diabat_projections[0,index] = energy
        diabat_projections[1,index] = diabat_projection
        #assiging the porjection and corresponding eigenstate energy to a column in the output array

    return diabat_projections

def trajectory_diabat_projection(trajectory_list, pseudoH_file_path, number_diabats, diabat_index):
    '''
    Function that essentially does what the above function does, but for a list of multiple trajectories

    Input: trajectory_list: iterable object filled with trajectory indices; pseudoH_file_path: string, path to pseudo-hamiltonian file; number_diabats: int, total number of diabats in active region; diabat_index: index of diabat whose projection you want (index starts at zero)

    Output: 2D array, where 1st row contains all eigenstates of the initial Hamiltonian of each trajectory included in the input list, second row contains all of these eigenstates' diabat projections. This means the list still has two rows, but the number of columns is now number_trajectories*number_diabats
    '''

    traj_diabat_projections = np.zeros((2, len(trajectory_list)*number_diabats ))
    counter = 0

    for trajectory in trajectory_list:

        diabat_projections = initial_diabat_projection(f'run-fssh-{trajectory}/{pseudoH_file_path}', number_diabats, diabat_index)
        traj_diabat_projections[:,counter: counter + number_diabats] = diabat_projections
        counter = counter + number_diabats

    return traj_diabat_projections


def IPR(population_frame, number_diabats, number_excitons):
    '''
    Function that calculates electron, hole and exciton IPR for a single snapshot of diabatic populations

    Input: population_frame: 1D array of diabats' populations wrt the surface/electronic wavefunction; number_diabats: int, number of diabats in the system; number_excitons: int, number of diabats that are excitonic

    Output: 3 floats corresponding to the IPR of the exciton, the acceptor carrier and the donor carrier
    '''

    exciton_populations = population_frame[-number_excitons:]
    exciton_populations = exciton_populations/np.sum(exciton_populations)
    #getting exciton diabatic populations and renormalising them

    exciton_IPR = 1/np.sum(exciton_populations**2)
    #then calculating exciton IPR from these normalised populations. This basically gives you the IPR if your wavefunction is only made up of excitons

    number_CT_states = (number_diabats - number_excitons)
    number_DA_molecules = int(np.sqrt(number_CT_states))
    #getting the number of donor/acceptor molecules from the number of CT-states. This assumes that there is an equal number of donors and acceptors in the system

    acceptor_site_pops = np.zeros(number_DA_molecules)
    donor_site_pops = np.zeros(number_DA_molecules)

    for j in range(number_DA_molecules):
    #i get the electron/hole IPR by defining new electronic states, which correspond to the hole localised on one molecule and the electron being wherever, and vice versa

        donor_sites = population_frame[j*number_DA_molecules : j*number_DA_molecules + number_DA_molecules]
        #get the expansion coeffs of all diabats where the donor's carrier is on the same donor molecule

        acceptor_sites = population_frame[np.arange(j, number_CT_states ,number_DA_molecules)]
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

    return exciton_IPR, acceptor_IPR, donor_IPR


def initial_IPR(sh_file_path, pseudoH_file_path, number_diabats, number_excitons):
    '''
    Function for calculating the exciton, hole and electron IPRs for the initial surface wavefunction of a single trajectory.

    Input: sh_file_path: string, path to run-sh-1.log file; pseudoH_file_path: string, path to run-pseudo-hamilt-1.xyz file; number_diabats: int, number of diabats composing the electronic wavefunction; number_excitons: int, out of the total number of diabats, how many of them are frenkel excitons? 

    Output: IPR values: floats, number of molecules over which the electron, hole and exciton is spread.
    '''

    sh_lines = xyz_parse_section(sh_file_path, 5, 1, 2)
    initial_adiabat = int(sh_lines[0][-1]) - 1
    #reading first line of the run-sh-1.log file, where it tells you your initial active state

    pseudoH_lines = xyz_parse_first_section(pseudoH_file_path, 3)
    H = build_sim_H(pseudoH_lines, number_diabats)
    #build the initial Hamiltonian

    eigenvals, eigenvecs = get_eigen(H)
    #then diagonalise to get all possible surface wavefunctions

    adiabat_wavefunction = eigenvecs[:,initial_adiabat]
    #the eigenvecs and active state labels are ordered with respect to energy, so you can use the active state label to index the correct eigenvector (i.e surface wavefunction)

    adiabat_populations = adiabat_wavefunction**2

    exciton_IPR, acceptor_IPR, donor_IPR = IPR(adiabat_populations, 420, 20)

    return exciton_IPR, acceptor_IPR, donor_IPR


def trajectory_initial_IPR(trajectory_list, sh_file_path, pseudoH_file_path, target_dir, number_diabats, number_excitons):
    '''
    Function for calculating the initial e/h/XT IPR of all trajectories specified in the list, and puts them in a 2D array.
    '''

    IPR_output_array = np.zeros((len(trajectory_list), 3))

    for index, trajectory in enumerate(trajectory_list):

        exciton_IPR, acceptor_IPR, donor_IPR = initial_IPR(f'{target_dir}/run-fssh-{trajectory}/{sh_file_path}', f'{target_dir}/run-fssh-{trajectory}/{pseudoH_file_path}', number_diabats, number_excitons)

        IPR_output_array[index, 0] = exciton_IPR
        IPR_output_array[index, 1] = acceptor_IPR
        IPR_output_array[index, 2] = donor_IPR

    return IPR_output_array

def single_trajectory_IPR(coeff_file_path, number_diabats, number_excitons):
    '''
    Function that calculated time-dependent IPR-vals for a single trajectory

    Input: coeff_file_path: str, path to trajectory's coeff file; number_diabats: int, number of diabats in the system; number_excitons: int, number of diabats that are excitons

    Output: 2D array of e/h/XT IPR-vals, 3 columns for the 3 types of IPR, number of rows = number of printed timesteps
    '''

    diabat_populations = diabatic_populations(coeff_file_path, number_diabats)
    IPR_array = np.zeros((len(diabat_populations), 3))

    for frame, population in enumerate(diabat_populations):

        exciton_IPR, acceptor_IPR, donor_IPR = IPR(population, number_diabats, number_excitons)

        IPR_array[frame, 0] = exciton_IPR
        IPR_array[frame, 1] = acceptor_IPR
        IPR_array[frame, 2] = donor_IPR

    return IPR_array

def diabat_locations(donor_coords, acceptor_coords, exciton_coords, population_frame, number_diabats, number_excitons):
    '''
    Function that calculates expectation values of particles' positions, for a single snapshot of diabatic populations

    Input: donor_coords/acceptor_coords: array(float), xyz-coms of all active donor/acceptor molecules; population_frame: array, diabat populations of the surface/electronic wavefunction; number_diabats: int, number of diabats in the system; number_excitons: int, number of diabats that are excitons

    Output: 3 floats that correspond to the exciton, acceptor carrier, and donor carrier expected positions
    '''

    exciton_populations = population_frame[-number_excitons:]
    exciton_populations = exciton_populations/np.sum(exciton_populations)
    #getting exciton diabatic populations and renormalising them

    weighted_exciton_locations = exciton_populations*exciton_coords
    #multiplying exciton diabat locations by the diabats' respective electronic populations

    average_exciton_location = np.sum(weighted_exciton_locations)
    #summing these weighted terms together to get the average position of the exciton

    number_CT_states = (number_diabats - number_excitons)
    number_DA_molecules = int(np.sqrt(number_CT_states))
    #getting the number of donor/acceptor molecules from the number of CT-states. This assumes that there is an equal number of donors and acceptors in the system

    acceptor_site_pops = np.zeros(number_DA_molecules)
    donor_site_pops = np.zeros(number_DA_molecules)

    for j in range(number_DA_molecules):
    #i get the electron/hole IPR by defining new electronic states, which correspond to the hole localised on one molecule and the electron being wherever, and vice versa

        donor_sites = population_frame[j*number_DA_molecules : j*number_DA_molecules + number_DA_molecules]
        #get the expansion coeffs of all diabats where the donor's carrier is on the same donor molecule

        acceptor_sites = population_frame[np.arange(j, number_CT_states ,number_DA_molecules)]
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

    weighted_donor_locations = donor_site_pops*donor_coords
    weighted_acceptor_locations = acceptor_site_pops*acceptor_coords

    average_donor_location = np.sum(weighted_donor_locations)
    average_acceptor_location = np.sum(weighted_acceptor_locations)
    #getting the average locations of the charge carriers in the donor and acceptor phase
   
    return average_exciton_location, average_acceptor_location, average_donor_location

 
def initial_diabat_locations(donor_coords, acceptor_coords, exciton_coords, sh_file_path, pseudoH_file_path, number_diabats, number_excitons):
    '''
    Function that uses the initial Hamiltonian and active surface of a trajectory to calculate the initial particle positions of that trajectory

    Input: donor_coords/acceptor_coords: array(float), xyz-coms of all active donor/acceptor molecules; sh_file_path: str, path to surface-hopping log file that contains the initial active surface; pseudoH_file_path: str, path to the pseudo-hamilt file that contains the initial Hamiltonian's elements; number_diabats: int, number of diabats in the system; number_excitons: int, number of diabats that are excitons

    Output: 3 floats that correspond to the exciton, acceptor carrier, and donor carrier expected positions
    '''

    sh_lines = xyz_parse_section(sh_file_path, 5, 1, 2)
    initial_adiabat = int(sh_lines[0][-1]) - 1
    #reading first line of the run-sh-1.log file, where it tells you your initial active state

    pseudoH_lines = xyz_parse_first_section(pseudoH_file_path, 3)
    H = build_sim_H(pseudoH_lines, number_diabats)
    #build the initial Hamiltonian

    eigenvals, eigenvecs = get_eigen(H)
    #then diagonalise to get all possible surface wavefunctions

    adiabat_wavefunction = eigenvecs[:,initial_adiabat]
    #the eigenvecs and active state labels are ordered with respect to energy, so you can use the active state label to index the correct eigenvector (i.e surface wavefunction)

    adiabat_populations = adiabat_wavefunction**2

    avg_XT_loc, avg_acc_loc, avg_don_loc = diabat_locations(donor_coords, acceptor_coords, exciton_coords, adiabat_populations, number_diabats, number_excitons)
   
    return avg_XT_loc, avg_acc_loc, avg_don_loc

def single_trajectory_diabat_locations(donor_coords, acceptor_coords, exciton_coords, coeff_file_path, number_diabats, number_excitons):
    '''
    Function that calculates time-dependent particle positions for a single trajectory

    Input: coeff_file_path: str, path to trajectory's coefficient file

    Output: 3D array of IPR-vals, each column is for a different type of expected location, whilst the number of rows = number of printed timesteps in the trajectory
    '''

    diabat_populations = diabatic_populations(coeff_file_path, number_diabats)
    #parse coeff file into array of diabatic populations for each timestep

    location_array = np.zeros((len(diabat_populations), 3))
    #initialise location array 

    for frame, population in enumerate(diabat_populations):
        #iterating over every timestep for which the diabatic populations were printed

            xt_loc, acc_loc, don_loc = diabat_locations(donor_coords, acceptor_coords, exciton_coords, population, number_diabats, number_excitons)
#using the diabatic populations for each timestep to get particles' expected locations for a single timestep

            location_array[frame, 0] = xt_loc
            location_array[frame, 1] = acc_loc
            location_array[frame, 2] = don_loc
            #then pass these into the array

    return location_array


