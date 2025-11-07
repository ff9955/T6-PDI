import numpy as np

def block_statistics(observable_array, number_blocks):
    
    block_means = np.zeros((number_blocks, len(observable_array[0])))

    number_trajectories = len(observable_array)
    trajectories_per_block = number_trajectories//number_blocks

    for block_index in range(number_blocks):
        block_means[block_index, :] = np.mean(observable_array[block_index*trajectories_per_block: block_index*trajectories_per_block + trajectories_per_block,:] ,axis = 0)

    final_block_mean = np.mean(block_means, axis = 0)
    final_block_sd = np.std(block_means, axis = 0)

    return final_block_mean, final_block_sd

def adjust_CS(CSS_single_populations, FS_single_populations, total_XT_single_populations, CT_single_populations, INT_single_populations, FS_pop_threshold):

    CSS_indices = np.where(FS_single_populations > FS_pop_threshold)

    for index in range(len(CSS_indices[0])):

        CSS_single_populations[CSS_indices[0][index], CSS_indices[1][index]:] = 1
        CT_single_populations[CSS_indices[0][index], CSS_indices[1][index]:] = 0
        INT_single_populations[CSS_indices[0][index], CSS_indices[1][index]:] = 0
        total_XT_single_populations[CSS_indices[0][index], CSS_indices[1][index]:] = 0

    return CSS_single_populations, total_XT_single_populations, CT_single_populations, INT_single_populations

def assign_full_separation(analysis_path, simulation, file_tag, FS_pop_threshold, save_path):

    total_XT_single_populations = np.loadtxt(f'{analysis_path}/{simulation}/analysis_arrays/{file_tag}_XT_individual_populations.txt')
    CSS_single_populations = np.loadtxt(f'{analysis_path}/{simulation}/analysis_arrays/{file_tag}_CSS_individual_populations.txt')
    INT_single_populations = np.loadtxt(f'{analysis_path}/{simulation}/analysis_arrays/{file_tag}_INT_individual_populations.txt')
    FS_single_populations = np.loadtxt(f'{analysis_path}/{simulation}/analysis_arrays/{file_tag}_FS_individual_populations.txt')
    CT_single_populations = 1 - INT_single_populations - total_XT_single_populations - CSS_single_populations

    CSS_single_populations, total_XT_single_populations, CT_single_populations, INT_single_populations = adjust_CS(CSS_single_populations, FS_single_populations, total_XT_single_populations, CT_single_populations, INT_single_populations, FS_pop_threshold)
    
    avg_CSS_pop, CSS_sd = block_statistics(CSS_single_populations, 5)
    avg_INT_pop, INT_sd = block_statistics(INT_single_populations, 5)
    avg_XT_pop, XT_sd = block_statistics(total_XT_single_populations, 5)
    avg_CT_pop, CT_sd = block_statistics(CT_single_populations, 5)

    states = ['XT', 'CT', 'INT', 'CSS']
    populations = [avg_XT_pop, avg_CT_pop, avg_INT_pop, avg_CSS_pop]
    for k in range(len(states)):

        np.savetxt(save_path + f'{simulation}_avg_{states[k]}.txt', populations[k])

def get_total_avg_property(ipr_array, population_array, pop_threshold):
    '''
    Function that computes total average of an observable using the electronic wavefunction, only when the population of a certain state-type
    (e.g XT/CSS) is above a threshold
    '''

    property_indices = np.where(population_array > pop_threshold)
    property_vals = ipr_array[property_indices]

    return np.mean(property_vals)

def band_stats(path, simulation, file_tag, state, stage):
    '''
    Function that returns mean and standard deviation of distribution of energies belonging to a particular eigenstate
    '''
    band_energies = np.loadtxt(f'{path}' + f'{simulation}/{file_tag}_{state}_{stage}DOS.txt')*27.211
    return np.mean(band_energies), np.std(band_energies)

def mechanistic_pathways(iCT_single_populations, total_XT_single_populations, total_CT_single_populations, CSS_single_populations, XT_CT_threshold, CSS_iCT_threshold, start, end):
    '''
    Function that tracks evolution of electronic wf over each printed timestep of every trajectory, categorises the current electronic state based on total diabatic populations of different 
    state types, and returns list of the pathway of states that each trajectory takes, expressed in string format

    Input: population arrays: 2D numpy array of floats, total diabatic populations of each state-type for different timesteps and trajectories;
    XT_CT_threshold, CSS_iCT_threshold: float, threshold of diabatic population for wf to be classified as a given state at a single timestep;
    start, end: int, specify the range of timesteps you want to sample in each trajectory

    Output: list of lists of strings, where each string element signposts the electronic state that domiinates the wf at that given timestep
    '''

    number_trajectories = len(CSS_single_populations)
    all_pathways = []
    #defining list that will contain the pathways of all trajectories

    for index in range(0,number_trajectories):

        pathway = []
        condensed_pathway = []

        for frame in range(start,end): #choosing timestep range to scan

            #for each printed timestep in the simulation, we categorise the wvaefunction into different types of electronic state,
            #depending on which type of diabat dominates the electronic wavefunction

            if CSS_single_populations[index][frame] > CSS_iCT_threshold:
                pathway.append('CSS')
                #e.g. if CSS diabats make up 80% of the wavefunction, then we label this timestep 'CSS', and we do this
                #for the other diabat types

            elif (iCT_single_populations[index][frame]) > CSS_iCT_threshold:
                pathway.append('iCT')
            elif total_CT_single_populations[index][frame] > XT_CT_threshold:
                pathway.append('CT')
            elif total_XT_single_populations[index][frame] > XT_CT_threshold:
                pathway.append('XT')

            else:

                pathway.append('H')
                #if neither CT not XT-states generally dominate the wavefunction, then we assign it as a hybrid state
            
            #   if total_XT_single_populations[index][frame] > population_threshold:
            #       pathway.append('hXT')

            #   elif total_CT_single_populations[index][frame] > population_threshold:
            #       pathway.append('hCT')

            #  else:
            #      pathway.append('H')

        condensed_pathway.append(pathway[0])
        for index2 in range(1, len(pathway)):
            #printing each timestep's wavefunction character is too complicated, so we just print the instances where they
            #change, to see how the carriers change throughout the simulation
        
            if pathway[index2] != pathway[index2 - 1]:
                condensed_pathway.append(pathway[index2])

        all_pathways.append(condensed_pathway)
        #then append this abridged pathway to the list of all trajectories' pathways

    return all_pathways

def group_pathways(all_pathways):
    '''
    Function that takes the list of lists describing the evolution of each trajectory's wf, gets rid of redundant sections, groups pathways 
    together of they're the same, and returns a list of printable pathways
    '''

    #fluctuation_counter = 0
    for index3 in range(len(all_pathways)):

        p = all_pathways[index3]
        fluctuation_keyword = False

        if (p.count('XT') > 1):
            #remove XT oscillating with other states before dissociating

            reversed_p = list(reversed(p))
            reversed_XT_index = reversed_p.index('XT')
            #find final occurence of XT

            final_XT_index = -reversed_XT_index - 1
            abridged_p = p[final_XT_index:]
            leftover_p = p[:final_XT_index]

            #if (leftover_p.count('iCT') > 0):
            #    fluctuation_counter += 1
            #    abridged_p = ['XT','H','iCT','H'] + abridged_p
            #    fluctuation_keyword = True

            #elif (leftover_p.count('CT') > 0):
            #    fluctuation_counter += 1
            #    abridged_p = ['XT','H','CT','H'] + abridged_p
            #    fluctuation_keyword = True

            p = abridged_p[:]

            #remove all elements before this occurence. Essentially, if a simulation begins in XT, we don't care about anything it
            #does before the exciton irreversibly dissociates

        if p[-1] == 'CSS':
            #remove CSS oscillations with other CT-state types

            CSS_index = p.index('CSS')
            #find first CSS instance

            p = p[:CSS_index+1] #+ [p[-1]]
            #remove all elements after this, we assume all charges are swept away by field once they reach the CSS

        if p[-1] == 'iCT':

            iCT_index = p.index('iCT') #find index of first iCT instance
            p = p[:iCT_index+1] #remove all elements after this to remove irrelevant fluctuations, and we know the traj. ends in the same state

        if p[-1] == 'CT':

            CT_index = p.index('CT')
            p = p[:CT_index+1] #removing all fluctuations after first non-interfacial CT occurrance

        all_pathways[index3] = p

    all_pathways = ['-->'.join(p) for p in all_pathways]
    #convert all truncated pathways to strings

    pathway_types = list(set(all_pathways))
    #define new list which contains strings of unique pathways only

    return all_pathways, pathway_types

def count_pathways_detailed(all_pathways, pathway_types):
    '''
    Function that counts number of trajectories that fall into a number of mechanistic categories, then returns the relative contributions 
    of these decay pathways

    Input: all_pathways: list of strings where each string illustrates the sequence of state-types that the electronic wf of one trajectory occupies 
    over time; pathway types: set formed from all_pathways list so it doesn't contain any duplicates

    Output: numpy array of floats, each element shows the relative contribution of a given decay pathway in comparison with all other pathways, except
    undissociated excitons since we're mostly interested in what happens after exciton dissociation.
    '''

    pathway_counts = map(lambda x: all_pathways.count(x), pathway_types)
    pathway_counts = list(pathway_counts)
    #get number of times these unique pathways occur in the trajectory ensemble

    pathway_1 = 0 #XT
    pathway_2 = 0 #XT -> H -> iCT
    pathway_3 = 0 #XT -> H -> niCT/CSS -> iCT
    pathway_4 = 0 #XT -> iCT -> niCT
    pathway_5 = 0 #XT -> H -> niCT

    #extra pathway for additional T6-PDI paper
    pathway_6 = 0 #XT -> H -> CSS -> niCT

    pathway_7 = 0 #XT -> H -> iCT -> CSS
    pathway_8 = 0 #XT -> H -> niCT -> CSS

    for index4 in range(len(pathway_types)):

        pathway_type = pathway_types[index4]
        pathway_count = pathway_counts[index4]

        #ignoring trajectories that aren't initialised as pure XT states, to stay consistent with decay pathways identified when plotting XT dissociation
        #distance distributions (natcomm project)
        #if pathway_type[0:2] != 'XT': continue

        #2nd T6-PDI paper: adding extra routes for XT dissociation to account for initial hybrid states instead of lumping them in with excitons, also adding 
        #routes to CSS straight from hybrid state and routes to niCT from CSS being populated first

        #ignoring initial XT and then hybrid states to see the relative contribution of each pathway for a given initial state
        #if pathway_type[0:2] != 'XT': continue

        if pathway_type[-3:] == 'CSS':

            if 'iCT' in pathway_type:
                pathway_7 += pathway_count
            else:
                pathway_8 += pathway_count
    
        elif pathway_type[-3:] == '>CT':

            if 'CSS' in pathway_type:
                pathway_6 += pathway_count
            elif 'iCT' in pathway_type:
                pathway_4 += pathway_count
            else:
                pathway_5 += pathway_count

        elif pathway_type[-3:] == 'iCT':
            if ('>CT' in pathway_type) or ('CSS' in pathway_type):
                pathway_3 += pathway_count
            else:
                pathway_2 += pathway_count

        else:
            pathway_1 += pathway_count

    #relative contributions of all pathways shown above
    all_types = np.array([pathway_1, pathway_2, pathway_3, pathway_4, pathway_5, pathway_6, pathway_7, pathway_8])/len(all_pathways)

    #only returning contributions of dissociated pathways, re-normalised to exclude undissociated excitons
    return all_types[1:]/np.sum(all_types[1:]), pathway_counts

def count_starting_points(all_pathways, pathway_types, pathway_counts):
    '''
    Function that counts the number of dissociated trajectories that originate from a given state - hybrid or full exciton
    '''
    #number of dissociated trajectories originating from hybrid/full XT states at t=0
    from_XT = 0
    from_H = 0
    still_XT = 0

    for index5 in range(len(pathway_types)):

        pathway_type = pathway_types[index5]
        pathway_count = pathway_counts[index5]

        if (pathway_type[-3:] != '>XT') and (pathway_type[-2:] != 'XT'):

            if pathway_type[0] == 'H':
                from_H += pathway_count
            elif pathway_type[0:2] == 'XT':
                from_XT += pathway_count
    
        else: still_XT += pathway_count

    return from_XT/len(all_pathways), from_H/len(all_pathways), still_XT/len(all_pathways)