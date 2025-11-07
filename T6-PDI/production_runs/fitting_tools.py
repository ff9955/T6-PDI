import numpy as np

def mexp_exciton_decay(time, initial_yval, parameter_list):

    decay_constant = parameter_list[0]

    return initial_yval*np.exp(-decay_constant*time)

def biexp_exciton_decay(time, initial_yval, parameter_list):
    #bi-exponential used to model exciton decay

    fast_decay, slow_decay = parameter_list[0], parameter_list[1]
    #decay constants associated with fast decay of near-interfacial excitons, and slower decay of non-interfacial excitons

    coeff_ratio = parameter_list[2]
    #ratio of slow to fast decay

    fast_coeff = initial_yval/(abs(coeff_ratio)+1) # x + y = i, y = rx, (1+r)x = i, x = i/(r+1), y = i - x
    slow_coeff = initial_yval - abs(fast_coeff)
    #getting the individual coeff values from the ratio and initial exciton population

    return fast_coeff*np.exp(-fast_decay*time) + slow_coeff*np.exp(-slow_decay*time)


def triexp_exciton_decay(time, initial_yval, parameter_list):
    #tri-expoential used to model exciton decay

    decay_1, decay_2, decay_3 = parameter_list[0], parameter_list[1], parameter_list[2]
    ratio_1, ratio_2 = parameter_list[3], parameter_list[4]

    coeff_1 = initial_yval/(1 + abs(ratio_1) + abs(ratio_2)) #x + y + z = i, y = rx, z=qx, x = i/(1+r+q)
    coeff_2 = initial_yval - abs(coeff_1) - abs(ratio_2*coeff_1) #y = i - x - z = i - x - qx
    coeff_3 = initial_yval - abs(coeff_1) - abs(coeff_2) #z = i - x - y

    return coeff_1*np.exp(-decay_1*time) + coeff_2*np.exp(-decay_2*time) + coeff_3*np.exp(-decay_3*time)


def residual(parameter_list, time, pop, model):
    #calculating linear residuals between data and model

    initial_pop = pop[0]
    return (pop - model(time, initial_pop, parameter_list))

def get_biexp_coeffs(ratio, initial_yval):

    fast_coeff = initial_yval/(ratio+1)
    slow_coeff = initial_yval - fast_coeff

    return fast_coeff, slow_coeff

def get_triexp_coeffs(ratio_1, ratio_2, initial_yval):

    coeff_1 = initial_yval/(1 + ratio_1 + ratio_2)
    coeff_2 = initial_yval - coeff_1 - ratio_2*coeff_1
    coeff_3 = initial_yval - coeff_1 - coeff_2

    return coeff_1, coeff_2, coeff_3