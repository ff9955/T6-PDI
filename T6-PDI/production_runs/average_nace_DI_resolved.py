import numpy as np

def accumulate_avg(data_array):

    number_elements = np.arange(1,len(data_array)+1)
    cumulative_sum = np.cumsum(data_array)

    return cumulative_sum/number_elements

initial_state = 'XT-H'
resolution_type = 'Ir'
resolution_dir = 'IPR_resolved'

file_tags = [f'{resolution_type}_{initial_state}_naces.txt', f'{resolution_type}_{initial_state}_weights.txt', f'{resolution_type}_{initial_state}_wnaces.txt', f'{resolution_type}_{initial_state}_weight_sums.txt', f'{resolution_type}_{initial_state}_mnaces.txt']

directory_list = ['physopt_retry', 'epsilon5_3xCT']
#obs_bins = np.arange(0,220,20)[1:]
obs_bins = np.arange(0,12,2)[1:]

resolved_avg_naces = np.zeros((len(directory_list), len(obs_bins)))
resolved_avg_weights = np.zeros((len(directory_list), len(obs_bins)))
resolved_avg_ediffs = np.zeros((len(directory_list), len(obs_bins)))
resolved_weight_sums = np.zeros((len(directory_list), len(obs_bins)))
resolved_avg_wnaces = np.zeros((len(directory_list), len(obs_bins)))
resolved_mnaces = np.zeros((len(directory_list), len(obs_bins)))

output_array_list = [resolved_avg_naces, resolved_avg_weights, resolved_avg_wnaces, resolved_weight_sums, resolved_avg_ediffs, resolved_mnaces]

for index in range(len(directory_list)):
    for index2 in range(len(file_tags)):

        for ob in range(len(obs_bins)):

            with open(f'../all_nace_files/{resolution_dir}/{initial_state}/{directory_list[index]}_naces/{directory_list[index]}_{obs_bins[ob]}_{file_tags[index2]}') as obs_file:

                obs_lines = obs_file.readlines()
                obs_filtered = filter(lambda x: x != 'None\n', obs_lines)
                obs_filtered = list(obs_filtered)
            
            formatted_obs = [obs[:-1] for obs in obs_filtered]
            observables = np.array(formatted_obs, dtype=float)
            print(len(observables))

#            acc_avg = accumulate_avg(abs(observables))
#            np.savetxt(f'{directory_list[index]}_{obs_bins[ob]}_acc_avg_{file_tags[index2]}', acc_avg)

            mean = np.mean(abs(observables))
            output_array_list[index2][index, ob] = mean

#np.savetxt(f'{initial_state}_avg_naces.txt', output_array_list[0])
#np.savetxt(f'{initial_state}_avg_weights.txt', output_array_list[1])
#np.savetxt(f'{initial_state}_avg_wnaces.txt', output_array_list[2])
#np.savetxt(f'{initial_state}_weight_sums.txt', output_array_list[3])
#np.savetxt(f'{initial_state}_ediffs.txt', output_array_list[4])
#np.savetxt(f'{initial_state}_max_naces.txt', output_array_list[-1])

print(output_array_list[0])
