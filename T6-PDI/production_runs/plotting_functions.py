import numpy as np
import matplotlib.pyplot as plt

def plot_population(ax_object, path, simulation, keyword_list, colour_list, label_list, pattern, thickness):

    for index in range(len(keyword_list)):

        population_array = np.loadtxt(path + f'/{simulation}_avg_{keyword_list[index]}.txt')
        ax_object.plot(np.arange(len(population_array))*10*0.05/1000, population_array, color=f'{colour_list[index]}', alpha=0.8, label=f'{label_list[index]}', linewidth=thickness, linestyle=pattern)

        if keyword_list[index] == 'XT':
            print(population_array[0])

def plot_hopping_model(ax_object, path, simulation, keyword_list, colour_list):

    for index in range(len(keyword_list)):

        population_array = np.loadtxt(path + f'/{simulation}_hm_{keyword_list[index]}_population.txt')
        ax_object.plot(np.arange(len(population_array))/100, population_array, color=f'{colour_list[index]}', alpha=0.8, linewidth=1.5, linestyle='--')

def exciton_decay(time, initial_yval, parameter_list):
    #bi-exponential used to model exciton decay

    fast_decay, slow_decay = parameter_list[0], parameter_list[1]
    #decay constants associated with fast decay of near-interfacial excitons, and slower decay of non-interfacial excitons

    coeff_ratio = parameter_list[2]
    #ratio of slow to fast decay

    fast_coeff = initial_yval/(coeff_ratio+1) # x + y = i, y = rx, (1+r)x = i, x = i/(r+1), y = i - x
    slow_coeff = initial_yval - fast_coeff
    #getting the individual coeff values from the ratio and initial exciton population

    return fast_coeff*np.exp(-fast_decay*time) + slow_coeff*np.exp(-slow_decay*time)

def plot_lsq_fit(ax_object, fitting_path, population_path, simulation, input_label=None):

    population_array = np.loadtxt(population_path + f'/{simulation}_avg_XT.txt')
    initial_population = population_array[0]
    simulation_time = np.arange(len(population_array))*10*0.05/1000

    fit_parameters = np.loadtxt(fitting_path + f'/{simulation}-fit.txt')
    ax_object.plot(simulation_time, exciton_decay(simulation_time, initial_population, fit_parameters[:4]), color = 'c', linestyle='--', label = input_label)

def plot_decay_constants(ax_object, fitting_path, simulation_list, parameter_index, colour_string, legend):

    parameter_list = []
    simulation_labels = ['(1)', '(2)', '(3)', '(4)', '(5)', '(6)', '(7)']
    
    for j in range(len(simulation_list)):

        fit_parameters = np.loadtxt(fitting_path + f'/{simulation_list[j]}-fit.txt')
        parameter_list.append(fit_parameters[parameter_index])

    ax_object.plot(simulation_labels, parameter_list, color=colour_string, linestyle='--', label=legend, alpha=0.9)
    ax_object.scatter(simulation_labels, parameter_list, color=colour_string, s=10, alpha=0.9)
    ax_object.set_yscale('log')

def plot_fitted_coefficients(ax_object, fitting_path, simulation_list, coeff_index, initial_yval_list):

    fast_coeffs = []
    slow_coeffs = []

    simulation_labels = ['(1)', '(2)', '(3)', '(4)', '(5)', '(6)', '(7)']
    
    for j in range(len(simulation_list)):

        fit_parameters = np.loadtxt(fitting_path + f'/{simulation_list[j]}-fit.txt')
        coeff_ratio = fit_parameters[coeff_index]
        initial_yval = initial_yval_list[j]

        fast_coeff = initial_yval/(coeff_ratio+1)
        slow_coeff = initial_yval - fast_coeff

        fast_coeffs.append(fast_coeff)
        slow_coeffs.append(slow_coeff)

    ax_object.plot(simulation_labels, slow_coeffs, color='purple', linestyle='--', alpha=0.9, label=r'$C_{\mathrm{b}}$')
    ax_object.scatter(simulation_labels, slow_coeffs, color='purple', s=10, alpha=0.9)

    ax_object.plot(simulation_labels, fast_coeffs, color='orange', linestyle='--', alpha=0.8, label=r'$C_{\mathrm{i}}$')
    ax_object.scatter(simulation_labels, fast_coeffs, color='orange', s=10, alpha=0.8)

def plot_DOS(ax_object, path, simulation, state_list, phase, colour_list, width_list, label_list, shading_list, bin_width_list):

    for index in range(len(state_list)):

        energy_array = np.loadtxt(path + f'{simulation}/{simulation}_{state_list[index]}_{phase}DOS.txt')*27.211
        ax_object.hist(energy_array, color=f'{colour_list[index]}', label=f'{label_list[index]}', rwidth=width_list[index], alpha=shading_list[index] ,orientation='horizontal', density=True,  bins=np.arange(-0.5,2.0,bin_width_list[index]))