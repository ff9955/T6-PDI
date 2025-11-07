import numpy as np
import matplotlib.pyplot as plt

def build_lattice(x_molecules, y_molecules, lattice_spacing):
    '''
    function builds a 1D/2D lattice with a prescribed number of evenly spaced molecules
    along the x/y-axes

    input: number of molecules along x and y-axes; lattice_spacing = distant (in angstroms)
    between adjacent molecules

    output: 2D array with two columns and n rows, each row contains x/y-coords
    of a single particle
    '''

    coordinate_list = []

    for index in range(x_molecules):
        for index2 in range(y_molecules):

            #nested loop that iterates over number of x and y-molecules
            #single x-mol taken, and iterated over all molecules over y-axes

            coordinate_list.append([index, index2])
            #appending list of x/y indeces of each particle, e.g [0,0] corresponds 
            #first y-molecule of first molecule along x-axis

    coordinate_vector = np.array(coordinate_list) #converting list of lists to 2D array
    coordinate_vector = coordinate_vector*lattice_spacing

    #distance between adjacent points is constant, so coordinates of each molecule are 
    #equal to indeces multiplied by spacing, e.g [0,1] = [0,5] as 0th x-molecule has
    #no distance, whereas 1st y-molecule is 5A away from 0th y-molecule
    

    return coordinate_vector


def DA_index(coord_vector, donors):
    '''
    function that sets indeces of donor/acceptor phase molecules using their coordinates

    input: 2D array of each molecules's coordinates, number of donors you want in your system

    output: 2 lists: one of the donor indeces and one with the acceptor indeces
    '''

    total_molecules = len(coord_vector)

    donor_indeces = [index for index in range(total_molecules)[:donors]]
    #donor indeces assinged by just taking the first 'donors' molecules in the lattice
    acceptor_indeces = [index2 for index2 in range(total_molecules)[donors:]]
    #acceptor indeces obtained by taking the remainder

    return donor_indeces, acceptor_indeces


def build_connectivity(coord_vector, max_distance):
    '''
    function that tells which molecules are nearest neighbours

    input: vector of lattice coordinates, distance cutoff beyond which points 
    are not nearest neighbours

    output: list of lists, where each list contains indeces of two molecules that 
    are closer than the distance cutoff
    '''

    connectivity_list = []
    
    for index in range(len(coord_vector)):
        #looping over all coordinates in coord_vector by indexing the array with 'index'

        if index == (len(coord_vector) - 1): break
        #breaks the for loop when you've reached the last element in range; avoids indexing error

        for index2 in range(index + 1, len(coord_vector)):
            #looping over the molecules that come after the molecule selected by 'index'

            euclidean_distance = np.sqrt((np.sum((coord_vector[index] - coord_vector[index2])**2)))
            #calculating distance between two molecules obtained from coord vector via 
            #'index' and 'index2' 

            if euclidean_distance <= max_distance:

                connectivity_list.append([index, index2])
                #if the two selected molecules are closer than the cutoff distance, a 2-component list
                #if appended to 'connectivity_list', where the components are the indeces of the
                #molecules in 'coord_vector'


    return connectivity_list

#it is important to understand that the indeces describing the molecules in the donor/acceptor
#lists and connectivity lists are the same


def define_CT_states_twophase(donor_list, acceptor_list):
    '''
    function that defines all possible charge transfer basis states of the system, operating under the assumption that electrons and holes
    are on distinct phases

    input: lists of indeces of molecules in donor and acceptor phases

    output: list of lists, where each list contains one donor and one acceptor index
    '''

    state_list = []

    for element in donor_list:
        for element2 in acceptor_list:
            state_list.append([element, element2])
    
    return state_list

def define_CT_states_single_phase(molecule_list):
    '''
    function that gets the CT-indices as the one above, but now under the assumption of a single phase, where the electrons and holes'
    locations are unrestricted

    input: list of all electronically active molecules

    output: list of lists, where each list contains index of molecule with the electron, and the index of the molecule with the hole
    '''

    state_list = []

    for hole_molecule in molecule_list:
        for electron_molecule in molecule_list:

            if hole_molecule != electron_molecule:
                state_list.append([hole_molecule, electron_molecule])
    
    return state_list

def define_XT_states(molecule_list):
    '''
    function that defines all possible exciton basis states of the system

    input: lists of indeces of molecules in donor or acceptor phase

    output: list of lists, where each list contains the index of a molecule twice,
    since both charge carriers are on the same molecule in a Frenkel exciton
    '''

    state_list = []

    for element in molecule_list:

        state_list.append([element, element])
    
    return state_list

def get_CT_distances(CT_states, coord_vector):
    '''
    Function that returns a list of x,y-distances between charge carriers of all CT-states in the system, in 
    the same order as the CT-state list passed into the function

    Input: List of CT-states (each CT-state is a 2-element list), array of coordinates of donor/acceptor points

    Output: lists of x-distances and y-distances between electrons and holes in each CT-state
    '''

    x_distance_list = []
    y_distance_list = []

    for CT in CT_states:

        donor = CT[0]
        acceptor = CT[1]

        xy_distances = coord_vector[acceptor] - coord_vector[donor]
        x_distance_list.append(xy_distances[0])
        y_distance_list.append(xy_distances[1])
    
    return x_distance_list, y_distance_list

def build_CT_block(state_list, coord_vector, connectivity, interaction_constant, donor_coupling, acceptor_coupling, e_potential):
    '''
    this function builds the CT block (site energies and couplings) of the electronic Hamiltonian

    input: list of CT-states, list of nearest neighbours, interaction constant of coulombic
    interaction between charge carriers, coupling magnitude between charge carriers, e_potential = electrochemical
    potential due to difference in charge carrier concentration, in units of V per angstrom

    output: NxN tight binding Hamiltonian matrix, where N = number of CT-states
    '''

    number_states = len(state_list)
    H_matrix = np.zeros((number_states, number_states))

    single_lattice_distance = np.sum(coord_vector[-1] - coord_vector[-2]) #the single-lattice distance is hard-coded and
    #depends here on the positions of the last 2 molecules in the array, remember this if you start to include non-uniform 
    #molecule separations

    e_potential = e_potential*1000*single_lattice_distance # conversion of electric field from V per angstrom
    # to meV per single lattice distance, no multiplication needed for conversion to eV because we're dealing
    # with single charge carriers

    for number in range(number_states):
        #looping over all states in state list

        for index in range(number_states)[number:]:
            #ignoring the states that come before 'index' in the states list

            if state_list[number] == state_list[index]:
                #if the states are the same, assign a site energy

                donor_molecule = state_list[number][0]
                acceptor_molecule = state_list[number][1]
                #obtain the indeces of the donor and acceptor molecules

                euclidean_distance = np.sqrt((np.sum((coord_vector[donor_molecule] - coord_vector[acceptor_molecule])**2)))
                #calculate distance between electron and hole

                euclidean_distance = euclidean_distance/single_lattice_distance
                #express this distance in units of lattice spacing

                H_matrix[number,index] = (interaction_constant)/euclidean_distance + e_potential*((-coord_vector[donor_molecule][0] + coord_vector[acceptor_molecule][0])/single_lattice_distance)
                #assign the site energy, which changes inversely with the distance

                H_matrix[number, index] = H_matrix[number, index]/2
                #divide by two as we're only assigning the top matrix triangle

            elif (0 in np.subtract(state_list[number], state_list[index])):
                #if the states are different, check that they are nearest neighbours before assigning an electronic coupling
                
                for element in connectivity:
                    if (state_list[number][0] == element[0]) and (state_list[index][0] == element[1]):
                        #this checks if the donor molecules are nearest neighbours
                        H_matrix[number, index] = donor_coupling
                        break

                    elif (state_list[number][1] == element[0]) and (state_list[index][1] == element[1]):
                        #this checks if the acceptor molecules are nearest neighbours
                        H_matrix[number, index] = acceptor_coupling
                        break
    
    H_matrix = H_matrix + H_matrix.T
    #assigning the bottom triangle with transpose of the top triangle, as the matrix is symmetric

    lowest_site_energy = min(np.diag(H_matrix))
    H_matrix = H_matrix - np.identity(len(H_matrix))*lowest_site_energy
    #subtracting all diagonal elements by site energy of interfacial CT-state, so all site energies are now
    #defined with the interfacial CT-state energy being the baseline

    return H_matrix


def build_XT_block(state_list, connectivity, offset, coupling_constant):
    '''
    this function builds the XT block (site energies and couplings) of the electronic Hamiltonian

    input: list of XT-states, list of nearest neighbours, energy offset of XT-states with respect to 
    interfacial CT-state (meV), coupling magnitude between charge carriers

    output: MxM tight binding Hamiltonian matrix, where M = number of XT-states
    '''

    number_states = len(state_list)
    H_matrix = np.zeros((number_states, number_states))
    
    for number in range(number_states):
        #looping over all states in state list

        for index in range(number_states)[number:]:
            #ignoring the states that come before 'index' in the states list

            if state_list[number] == state_list[index]: 
               H_matrix[number, index] = offset/2
            #if the states are the same, assign a site energy
            
            else:
                for element in connectivity:
                    
                    # if the states are different, check that they are nearest neighbours before assigning an electronic coupling

                    if (state_list[number][0] in element) and (state_list[index][0] in element):

                        H_matrix[number, index] = coupling_constant
                        break

    H_matrix = H_matrix + H_matrix.T

    return H_matrix


def build_XT_CT_block(XT_states, CT_states, connectivity, coupling=None, XT_hole_coupling=None, XT_electron_coupling=None):
    '''
    function that assigns the offdiagonal block of the full electronic Hamiltonian
    matrix, which describes coupling between XT and CT states, generalised to work on two-phase system with exciton on 1 phase (set 2nd/3rd
    coupling args to none) or single-phase where exciton is on every state and can couple to holes/electrons (assign values to the 2nd/3rd 
    coupling arguments)

    input: list of CT/XT states, list of nearest neighbours, coupling constant

    output: NxM matrix, where N = number of CT-states, and M = number of XT states
    '''

    if coupling==None:
        if XT_hole_coupling==None or XT_electron_coupling==None:
            raise ValueError('You must pass a coupling value into the function')
    else:
        if (XT_electron_coupling) or (XT_hole_coupling):
            raise ValueError('You must choose between a single or two-phase system')

    number_XT_states = len(XT_states)
    number_CT_states = len(CT_states)

    H_matrix = np.zeros((number_CT_states, number_XT_states))

    #if (XT_states[0][0] == CT_states[0][0]):
        #this checks if the excitons are in the donor or acceptor phase
    #    CT_index = 1
    #else:
    #    CT_index = 0

    #reminder: donor molecules are first in the state list, acceptors are second

    for number in range(number_CT_states):

        for index in range(number_XT_states):

            index_diff = np.subtract(CT_states[number], XT_states[index])

            if (0 in index_diff):
                #this checks that the non-nearest neighbour indeces are equal, since
                # coupling between charge carriers on four different sites is zero
                matching_loc = np.where(index_diff == 0)[0][0]

                for element in connectivity:
                    if (CT_states[number][matching_loc-1] in element) and (XT_states[index][matching_loc-1] in element):
                        #if either electrons or holes are on nearest neighbour sites, coupling is assigned

                        if (XT_electron_coupling!=None) and (XT_hole_coupling!=None):
                            #activates if you've passed in separate values for XT-e and XT-h coupling, which corresponds to single-phase system

                            if matching_loc == 0:
                                #matching loc == 0 means exciton on same site as hole and couples to nearest-neighbour electron
                                H_matrix[number, index] = XT_electron_coupling
                            else:
                                H_matrix[number, index] = XT_hole_coupling
                        else:
                            #just one coupling passed in --> two-phase system where exciton is on 1 phase only
                            H_matrix[number, index] = coupling
    
    return H_matrix


def build_full_Hamiltonian(CT_block, XT_block, XT_CT_block):
    '''
    function that combines CT, XT, and CT-XT blocks to create full state-space 
    Hamiltonian matrix where the CT and XT blocks are ordered in the same way each time

    input: different matrix blocks as defined above

    output: DxD matrix where D = total number of states
    '''

    number_CT_states = len(CT_block)
    number_XT_states = len(XT_block)
    total_states = number_CT_states + number_XT_states

    H_matrix = np.zeros((total_states, total_states))

    H_matrix[0:number_CT_states, 0:number_CT_states] = CT_block
    H_matrix[number_CT_states:total_states, number_CT_states:total_states] = XT_block
    #assigning the two diagonal blocks
    
    H_matrix[0:number_CT_states, number_CT_states:total_states] = XT_CT_block
    #assign the top-right XT-CT block

    H_matrix[number_CT_states:total_states, 0:number_CT_states] = XT_CT_block.T
    #assign the bottom-left XT-CT block by taking the transpose of the one above

    return H_matrix


def get_eigen(matrix):
    '''
    function computes energies and eigenvectors of eigenstates in increasing order of energy

    input: square Hamiltonian matrix

    output: row vector (length D) of eigen energies, matrix (DxD) of orthonormal eigenvectors
    where each column is a single eigenvector
    '''
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    idx = eigenvalues.argsort()
    #print idx
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    return eigenvalues, eigenvectors

def CT_site_populations_2phase(eigenvector, donors, acceptors):
    '''
    function that returns the distribution of a single charge carrier's eigenfunction
    across the system's lattice sites

    input: single eigenvector of Hamiltonian, lists of donor/acceptor indeces

    output: list of length equal to number of lattice sites; each element corresponds to
    the proportion of the eigenfunction that is present on this lattice site
    '''

    number_CT_states = len(donors)*len(acceptors)
    CT_states_per_donor = int(number_CT_states/len(donors))
    #each donor contributes the same number of CT-states to the overall number of CT-states

    CT_eigenvector = eigenvector[:number_CT_states]
    #remove the elements of the total eigenvector that correspond to the expansion coefficients
    #of the excitonic basis states
    donor_site_populations = []

    for index in range(0, number_CT_states, CT_states_per_donor):
        #iterate over list of CT-states, taking slices at a time, where each slice contains all
        #possible CT-states of a single donor molecule

        cumulative_population = np.sum(CT_eigenvector[index: index + CT_states_per_donor]**2)
        #we sum the populations of these CT-states, this tells us the total population that is present
        #on this one donor molecule
        donor_site_populations.append(cumulative_population)

    acceptor_site_populations = []
    cumulative_population = []
    #same principle holds for acceptor, although we index differently, as the state lists
    #were made by summing over all acceptor molecules per donor

    for acceptor_number in range(0, len(acceptors)):

        single_acceptor_populations = [CT_eigenvector[number2]**2 for number2 in range(0+acceptor_number, number_CT_states, CT_states_per_donor)]
        #list comprehension gets all CT-states corresponding to a single donor by skipping N elements at a time, where N = number of CT-states per donor
        #each donor has a CT-state with just one acceptor, so an acceptor's adjacent CT-states are always separated by N list elements

        cumulative_population = np.sum(single_acceptor_populations)
        #sum populations of all CT-states of a single acceptor

        acceptor_site_populations.append(cumulative_population)

    return donor_site_populations, acceptor_site_populations 

def CT_site_populations_1phase(eigenvector, CT_states, molecule_indices):

    hole_sites = np.array([element[0] for element in CT_states])
    electron_sites = np.array([element[1] for element in CT_states])

    hole_populations = []
    electron_populations = []

    for index in molecule_indices:

        hole_indices = np.where(hole_sites == index)[0]
        electron_indices = np.where(electron_sites == index)[0]

        cumulative_h_population = np.sum(eigenvector[hole_indices]**2)
        hole_populations.append(cumulative_h_population)

        cumulative_e_population = np.sum(eigenvector[electron_indices]**2)
        electron_populations.append(cumulative_e_population)
    
    return hole_populations, electron_populations


def XT_site_populations(eigenvector, XT_states):
    '''
    function which calculates how the excitonic populations of a charge carrier's eigenfunction is distributed over 
    a lattice

    input: eigenvector of a Hamiltonian, list of excitonic basis states

    output: list where each element shows the proportion of the eigenfunction that is present on a lattice
    site as an exciton
    '''

    number_XT_states = len(XT_states)
    XT_eigenvector = eigenvector[-number_XT_states:]

    site_populations = []

    for coefficient in XT_eigenvector:
        site_populations.append(coefficient**2)
    
    #simply summing over the squared elements of the eigenvector's XT-section here, since
    #each lattice site has just one possible exciton

    return site_populations


def integrated_CT_populations(donor_site_populations, acceptor_site_populations, y_molecules):
    '''
    function that sums the lattice sites' relative populations along the y-axis to show how
    an eigenfunction's CT-component is distributed over columns of the lattice, not only
    individual sites

    input: list of CT-state populations of donor and acceptor sites; number of molecules per 
    column in the lattice (y-axis)

    output: lists with length equal to half the number of molecules along the x-axis, each
    element corresponds to the proportion of the eigenfunction's CT-component on a certain column 
    '''

    integrated_donor_populations = [np.sum(donor_site_populations[number:number+y_molecules]) for number in range(0,len(donor_site_populations),y_molecules)]

    #these lists comprehensions just iterate over the per-site population lists, and sum over every N molecules, where N = the number of molecules in a single column

    integrated_acceptor_populations = [np.sum(acceptor_site_populations[number2:number2+y_molecules]) for number2 in range(0,len(acceptor_site_populations),y_molecules)]

    return integrated_donor_populations, integrated_acceptor_populations


def integrated_XT_populations(XT_populations, y_molecules):
    '''
    function that sums over relative XT-populations for each column of the lattice - same as above

    input: list of relative presence of eigenfunction's XT-component per site, number of molecules
    in a single column of the lattice

    output: lists with length equal to half the number of molecules along the x-axis, each
    element corresponds to the proportion of the eigenfunction's XT-component on a certain column 
    '''

    integrated_populations = [np.sum([XT_populations[number:number+y_molecules]]) for number in range(0, len(XT_populations), y_molecules)]

    return integrated_populations

def bin_eigenspectrum(eigenvalues, eigenvectors, energy_spacing):
    '''
    function that will sort the eigenvalues and eigenvectors obtained from an electronic Hamiltonian
    into different bins delineated by values of the eigen energies.

    Input: 1D array of eigenvalues, matrix of corresponding eigenvectors, the energy separation you want to leave
    between adjacent bins: the smaller the separation, the more bins you will have

    Output: separate lists and matrices of eigenvectors and eigenvalues, that fall within a certain bin of
    energy
    '''

    energy_boundaries = np.arange(eigenvalues[0], eigenvalues[-1] + energy_spacing, energy_spacing)
    #using the energy_spacing to assign the upper/lower bounds of each bin

    binned_eigenvalues = []
    binned_eigenvectors = []
    used_bins = []

    for index in range(len(energy_boundaries)-1):

        binned_eigenvals = filter(lambda x: energy_boundaries[index] <= x < energy_boundaries[index+1], eigenvalues)
        binned_eigenvals = list(binned_eigenvals)
        #using the filter function to obtain all eigenvalues that fall within a single bin's bounds

        binned_indices = [np.where(eigenvalues==val)[0][0] for val in binned_eigenvals]
        #getting the corresponding indices of these binned eigenvalues in the total eigenvalue list

        if not binned_eigenvals:
            continue

        used_bins.append([energy_boundaries[index], energy_boundaries[index+1]])

        number_binned_vals = len(binned_indices)
        eigenval_array = np.zeros(number_binned_vals)
        eigenvec_array = np.zeros((len(eigenvalues), number_binned_vals))

        for index2, binned_index in enumerate(binned_indices):

            eigenval_array[index2] = eigenvalues[binned_index]
            #making an array of binned eigenvals by indexing the chosen eigenvals from the full list

            eigenvec_array[:,index2] = eigenvectors[:, binned_index]
            #then doing the same with different columns of the eigenvector matrix, since the order of 
            #columns(eigenvectors) is the same as the order of eigenvalues
        
        binned_eigenvalues.append(eigenval_array)
        binned_eigenvectors.append(eigenvec_array)
        #appending one bin to total lists of arrays of binned eigenvalues and eigenvectors

    return binned_eigenvalues, binned_eigenvectors, used_bins

def average_binned_populations(binned_eigenvecs, donor_list, acceptor_list, exciton_states, y_molecules):
    '''
    function that will take the binned matrices of eigenvectors, calculate the XT and CT populations of each
    eigenvector, then average the populations for each bin of eigenvectors
    '''

    number_bins = len(binned_eigenvecs)
    x_donors = len(donor_list)//y_molecules
    x_acceptors = len(acceptor_list)//y_molecules
    x_excitons = len(exciton_states)//y_molecules

    donor_avg_population_array = np.zeros((number_bins, x_donors))
    acceptor_avg_population_array = np.zeros((number_bins, x_acceptors))
    exciton_avg_population_array = np.zeros((number_bins, x_excitons))

    for index in range(len(binned_eigenvecs)):

        matrix_transpose = binned_eigenvecs[index].T

        donor_population_array = np.zeros((len(matrix_transpose), x_donors))
        acceptor_population_array = np.zeros((len(matrix_transpose), x_acceptors))
        exciton_population_array = np.zeros((len(matrix_transpose), x_excitons))

        for index2 in range(len(matrix_transpose)):

            donor_populations, acceptor_populations = CT_site_populations_2phase(matrix_transpose[index2], donor_list, acceptor_list)
            exciton_populations = XT_site_populations(matrix_transpose[index2], exciton_states)

            integrated_donors, integrated_acceptors = integrated_CT_populations(donor_populations, acceptor_populations, y_molecules)
            integrated_excitons = integrated_XT_populations(exciton_populations, y_molecules)

            donor_population_array[index2,:] = integrated_donors
            acceptor_population_array[index2,:] = integrated_acceptors
            exciton_population_array[index2,:] = integrated_excitons
        
        donor_population_average = np.mean(donor_population_array, axis=0)
        acceptor_population_average = np.mean(acceptor_population_array, axis=0)
        exciton_population_average = np.mean(exciton_population_array, axis=0)

        donor_avg_population_array[index,:] = donor_population_average
        acceptor_avg_population_array[index,:] = acceptor_population_average
        exciton_avg_population_array[index,:] = exciton_population_average

    return donor_avg_population_array, acceptor_avg_population_array, exciton_avg_population_array

def plot_eigenspectrum(donor_populations, acceptor_populations, exciton_populations, figure_size, max_shading, exciton_radius=50, bins=None, eigenvals=None, active_state_energy=None, binned_eigenvals=None):
    '''
    function for plotting the relative diabat populations on a coloured grid, using either averaged populations of populations of single 
    eigenstates

    Input: donor/acceptor/exciton_populations: 2D arrays where the rows refer to diffrent molecules' population, and the columns refer to
    different eigenstates, max_shading refers to maximum level of shading on each square, this should be smaller if there is less variation
    in colour between squares, exciton_radius controls the size of the dots representing the exciton populations, bins is a list of lists,
    where each list has the bounds of binned eigenstates, whilst eigenvals is just a list of eigenvalues corresponding to the eigenvectors

    Output: 2D plot of coloured rows, where shading shows relative populations of holes/electrons on each site, and dot size shows the same
    thing but for the excitons
    '''

    number_eigenfunctions = len(donor_populations)
    fig, axs = plt.subplots(nrows = number_eigenfunctions, ncols = 2, sharex = True, figsize = figure_size)

    for band in range(number_eigenfunctions):
        #iterating over the number of distinct populations you're going to plot

        axs[band,0].imshow(np.array([donor_populations[-(band+1)]]), cmap = 'Reds', vmin = 0, vmax = max_shading)
        axs[band,0].xaxis.set_visible(False)
        axs[band,0].set_yticks([0]) #using the 'band' index to get the highest energy values first, because the graph is plotted from
        #top to bottom, and we want the lowest-energy states at the bottom for clarity.

        if bins:
            axs[band,0].set_yticklabels([str((round(bins[-(band+1)][0],1), round(bins[-(band+1)][-1],1)))])
        elif not type(eigenvals) == None:
            axs[band,0].set_yticklabels([str(round(eigenvals[-(band+1)],1))])
        else:
            raise ValueError('No eigenvalues or bins provided')
        #setting labels of each eigenstate population plot to be either the bounds of the bin or the eigen energy, if you don't pass either
        #one in, the function will throw an error

        axs[band,1].imshow(np.array([acceptor_populations[-(band+1)]]), cmap = 'Blues', vmin = 0, vmax = max_shading)
        axs[band,1].yaxis.set_visible(False)
        axs[band,1].xaxis.set_visible(False)
        #doing same thing but for the acceptor phase

        if binned_eigenvals: #if using binned eigenstates, adding the number of eigenstates used for each subplot
            number_binned_eigenstates = len(binned_eigenvals[-(band+1)])
            axs[band,1].text(1.03, 0.5, f'{number_binned_eigenstates}', transform = axs[band,1].transAxes, va = 'center')

        if active_state_energy:
            if bins:
                if ( bins[-(band+1)][0] < active_state_energy < bins[-(band+1)][-1]):
                    axs[band,1].axhline(0, color = 'green')
                elif active_state_energy == bins[-(band+1)][0]:
                    axs[band,1].axhline(0, color = 'green')
            elif eigenvals:
                if (active_state_energy == eigenvals[-(band+1)]):
                    axs[band,1].axhline(0, color = 'green')
        #putting a green streak through the subplot that contains the active eigenstate (if there is one)

        for i in range(len(exciton_populations[0])):

            axs[band,1].scatter(i, 0, c = 'm', s = int(exciton_radius*exciton_populations[-(band+1)][i]))
            #plotting the exciton populations as dots on the acceptor phase, this should also be generalised at some point

    fig.supxlabel('Donor                            Acceptor')

    return fig, axs

def plot_eigenspectrum_1phase(hole_populations, electron_populations, exciton_populations, eigenvals, figure_size, default_radius=300):

    number_eigenvectors = len(hole_populations)
    number_molecules = len(exciton_populations[0])
    
    fig, axs = plt.subplots(nrows=number_eigenvectors, ncols=1, figsize=figure_size)

    for band in range(number_eigenvectors):
        axs[band].set_yticks([0.5])
        axs[band].set_xticks([])
        axs[band].set_yticklabels([str(round(eigenvals[-(band+1)],1))])

        for molecule in range(number_molecules):

            axs[band].scatter(molecule, 0.5, s = default_radius*exciton_populations[-(band+1), molecule], color='k', alpha=0.80)
            axs[band].scatter(molecule, 0.5, s = default_radius*hole_populations[-(band+1), molecule], color='c', alpha=0.65)
            axs[band].scatter(molecule, 0.5, s = default_radius*electron_populations[-(band+1), molecule], alpha=0.70, facecolors='none', edgecolors='r')

    fig.supylabel('Eigenstate energy /meV')
    fig.tight_layout()