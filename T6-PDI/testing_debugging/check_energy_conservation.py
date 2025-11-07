#!/bin/python

# IMPORT
import re
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import os


windows_path = '/mnt/c/Users/samue/Desktop/Research_work/'

datadict = {
    
    "AOM_ANT_NACV_RESCAL" : "./AOM_ANT_NACV_RESCAL/data-general-AOM_ANT_NACV_RESCAL-2002210840",
}

number_atoms_mol = 24
number_active_mol = 36 
keys = ('TEMPERATURE',)


temperatures = {}
energy_drift = {}
for directory in datadict:
    dataname = datadict[directory]
    temperatures[directory] = {}
    energy_drift[directory] = {}
    for property in ['Temperature', 'Total-energy']:
        filename =  dataname + '/' + 'Mean-' + property + '-' + '-'.join(keys) + '.dat'
        prop = []
        for line in open(filename).readlines()[1:]:
            prop.append( line.split() )
        prop = np.array(prop).transpose() 
        if property == 'Temperature':
            temperatures[directory][keys] = np.mean( [ float(x) for x in prop[1] ] )
        elif property == 'Total-energy':
            energy_drift[directory][keys]  = np.mean( [abs(float(x)) for x in prop[3]])

totalener = []
temp_sim = []
temp_target = sort(datadict.keys())
print "| Energy drift (Ha/atom/ps) | Temperature (K)"
for temp in temp_target:
    totalener.append(energy_drift[temp][keys]  * 1000 / (number_active_mol*number_atoms_mol))
    temp_sim.append(temperatures[temp][keys])
    print "| %s  | %.1E | %.1F " % (temp, totalener[-1], temp_sim[-1])
