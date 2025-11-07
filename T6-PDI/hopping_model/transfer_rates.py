import numpy as np
from scipy.linalg import expm

def marcus_sc(reorgE, deltaA, V, freq, lowerE, upperE):
    '''
    Function for calculating semi-classical Marcus rate of charge transfer between diabats - see Jochen's 2015 review.
    Applies to both adiabatic and non-adiabatic limits of charge transfer. Estimating activation barriers by calculating 
    free energies of ground/excited state adiabats of 2-state system.

    Input: reorgE: float, reorganisation energy of transfer (in meV), deltaA: float, difference between free energy minima 
    of diabatic states (meV), should be passed in as Initial_state - Final_state, V: float, electronic coupling between states (meV).
    freq: float, frequency of effective nuclear mode (in wavenumbers); lowerE/upperE: float, upper and lower energy bounds for plotting
    the potential energy surfaces of diabats/adiabats

    Output: k_c: float, rate constant of charge transfer
    '''

    freq = freq*100*3e8 #conversion to Hz from wavenumbers
    kB = 1.38e-23*6.2415e18*1000 #meV K-1
    planck = 6.626e-34*6.2415e18*1000 #meV s-1

    deltaE = np.arange(lowerE, upperE)

    #free energies of diabats 1 and 2
    Aa = ((deltaE - (deltaA + reorgE))**2)/(4*reorgE)
    Ab = ((deltaE - (deltaA - reorgE))**2)/(4*reorgE) + deltaA

    #free energies of ground and excited state adiabats
    A0 = 0.5*(Aa + Ab) - 0.5*np.sqrt(deltaE**2 + 4*V)
    A1 = 0.5*(Aa + Ab) + 0.5*np.sqrt(deltaE**2 + 4*V)

    #index of transition state energy always where deltaE = 0
    ts_index = np.where(deltaE == 0)[0][0]

    #inverted marcus regime, hopping downwards in energy
    if deltaA <= -reorgE:
        #higher energy diabat starts in bottom of A1 P.E.S
        A1_min = np.min(A1)
        ts_energy = A1[ts_index] #T.S where diabatic energies intersect

        deltaA_dagger = max([ts_energy - A1_min, 0])
    
    #reverse reaction of the process above
    elif deltaA >= reorgE:

        A0_min = np.min(A0)
        ts_energy = A0[ts_index]

        #if transfer from A1 has no barrier, then this will approximately be equal to driving force
        deltaA_dagger = ts_energy - A0_min

    #normal Marcus regime, only ground state adiabat is now relevant
    #driving force can be positive or negative, we plot full potential energy surface (PES) and have to determine which minimum to take for the given transfer
    #reaction in the rate matrix we have to calculate

    #upward transition in energy, take RHS of PES
    elif deltaA > 0:
        
        ts_energy = A0[ts_index]
        A0_left = A0[:ts_index]
        A0_right = A0[ts_index+1:]
        A0_min = np.min(A0_right)

        #first check LHS to see if there's actually an activation barrier
        if np.min(A0_left) < ts_energy: #if barrier, use ts energy to get peak
            deltaA_dagger = ts_energy - A0_min
        else: #no barrier, use driving force
            deltaA_dagger = deltaA

    #downward transiton --> use RHS as well, but getting activation E is different
    elif deltaA < 0:

        ts_energy = A0[ts_index]
        A0_right = A0[ts_index+1:]
        A0_min = np.min(A0_right)

        #take energy difference between min and ts energy, or zero if 0/negative
        deltaA_dagger = max([ts_energy - A0_min, 0])

    #symmetric PES, only one side of barrier is sufficient to get activation E; putting in extra condition for clarity
    elif deltaA == 0:

        ts_energy = A0[ts_index]
        A0_right = A0[ts_index:]
        A0_min = np.min(A0_right)

        deltaA_dagger = max([ts_energy - A0_min, 0])

    else:
        raise ValueError('Invalid driving force specified')
    
    lz_factor = (V*(np.pi**1.5))/(planck*freq*np.sqrt(300*reorgE*kB))
    #2 pi gamma factor

    P_lz = 1 - np.exp(-lz_factor)
    #transition probability

    if deltaA < -reorgE:
        k_el = 2*P_lz*(1-P_lz)

    else:
        k_el = 2*P_lz/(1+P_lz)
    #transmission coefficient depends on relative values of reorgE and lz factor

    k_c = k_el*freq*np.exp(-deltaA_dagger/(300*kB))
    #put terms into total expression

    return k_c, deltaA_dagger


def marcus_sc_warshel(reorgE, deltaA, V, freq):
    '''
    Function for calculating semi-classical Marcus rate of charge transfer between diabats - see Jochen's 2015 review.
    Applies to both adiabatic and non-adiabatic limits of charge transfer. Using Warshel's equation for free energy of 
    activation (divergent for negative driving forces close in magnitude to reorganisation energy).

    Input: reorgE: float, reorganisation energy of transfer (in meV), deltaA: float, difference between free energy minima 
    of diabatic states (meV), should be passed in as Initial_state - Final_state, V: float, electronic coupling between states (meV).
    freq: float, frequency of effective nuclear mode (in wavenumbers)

    Output: k_c: float, rate constant of charge transfer
    '''

    freq = freq*100*3e8 #conversion to Hz from wavenumbers
    kB = 1.38e-23*6.2415e18*1000 #meV K-1
    planck = 6.626e-34*6.2415e18*1000 #meV s-1

    if abs(deltaA) < reorgE:
        deltaA_dagger = ((reorgE + deltaA)**2)/(4*reorgE) - np.sqrt(V) + V/(deltaA + reorgE)
    #inverted marcus regime
    elif (deltaA < -reorgE):
        deltaA_dagger = 0.0
    #driving force is positive and > reorgE --> dominates activation barrier
    elif (deltaA > 0) and (np.sqrt(V) >= reorgE/2):
        deltaA_dagger = deltaA
    elif (deltaA < 0) and (np.sqrt(V) >= reorgE/2):
        deltaA_dagger = 0
    else:
        deltaA_dagger = deltaA
    #getting activation free energy from deltaA

    lz_factor = (V*(np.pi**1.5))/(planck*freq*np.sqrt(300*reorgE*kB))
    #2 pi gamma factor

    P_lz = 1 - np.exp(-lz_factor)
    #transition probability

    if deltaA < -reorgE:
        k_el = 2*P_lz*(1-P_lz)

    else:
        k_el = 2*P_lz/(1+P_lz)
    #transmission coefficient depends on relative values of reorgE and lz factor

    k_c = k_el*freq*np.exp(-deltaA_dagger/(300*kB))
    #put terms into total expression

    return k_c, deltaA_dagger


def marcus_na(reorgE, deltaA, V):
    #marcus rate equation in nonadiabatic (small coupling) limit
    
    kB = 1.38e-23*6.2415e18*1000 #meV K-1
    planck = 6.626e-34*6.2415e18*1000 #meV s-1

    deltaA_dagger_na = ((reorgE + deltaA)**2)/(4*reorgE)
    fc = (1/np.sqrt(4*np.pi*reorgE*kB*300))*np.exp(-deltaA_dagger_na/(300*kB))

    k_c = (4*np.pi*np.pi/planck)*V*fc

    return k_c


def time_dependent_populations(rate_matrix, times, p_init):
    '''
    Function for getting the populations of all states in the rate matrix at all times passed into the function.

    Input: rate_matrix: 2D np array of floats, rate constants of transfer processes between different sites. times: 1D np array of floats,
    array of times for which one wants to calculate site populations. p_init: 1D float array, site populations at t=0.

    Output: p_t: 2D float matrix, each column corresponds to the populations of all sites at a single timestep.
    '''

    p_t = np.zeros((len(p_init), len(times)))
    p_t[:,[0]] = p_init
    #1st row of p_t will be the initial populations (t=0)

    #iterating through the non-zero timesteps
    for index in range(1,len(times)):
        t = times[index]

        #multiplying rate matrix by t
        kt = rate_matrix*t
        #exponentiating rate matrix
        kt = expm(kt)

        #populations for this time, t, obtained by multiplying exp(kt) matrix by p_init vector
        p_t[:,[index]] = np.matmul(kt,p_init)

    return p_t