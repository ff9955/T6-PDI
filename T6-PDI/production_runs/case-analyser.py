import os
import shutil
import time
import numpy as np
import sys

path_nve = sys.argv[1]

####input parameters#####
tot_traj = 100
#########################

# copy input files and create FSSH-* directories 
#for i in range(tot_traj):
#    path = path_nve + 'FSSH-' + str(i+1)
#    if not (os.path.isdir(path)):
#        os.mkdir(path)
#pop_ = []
pop1_ = 0
pop2_ = 0
pop3_ = 0
pop4_ = 0
pop5_ = 0
pop6_ = 0
for i in range(tot_traj):
    if os.path.isdir(f"{path_nve}/run-fssh-" + str(i)):
        path = path_nve + 'run-fssh-' + str(i) + '/' + 'cases.log'
        with open (path, 'r') as traj:
            all_lines = traj.readlines()
            for line in all_lines:
                if(line.find('CASE') != -1):
                    #pop_.append(int(line[line.find('CASE')+4]))
                    x = int(line[line.find('CASE')+4])
                    if (x == 1):
                        pop1_ = pop1_ + 1
                    elif (x == 2):
                        pop2_ = pop2_ + 1
                    elif (x == 3):
                        pop3_ = pop3_ + 1
                    elif (x == 4):
                        pop4_ = pop4_ + 1
                    elif (x == 5):
                        pop5_ = pop5_ + 1
                    elif (x == 6):
                        pop6_ = pop6_ + 1

t_pop = pop1_+pop2_+pop3_+pop4_+pop5_+pop6_

rate1 = float(pop1_)/float(t_pop)
rate2 = float(pop2_)/float(t_pop)
rate3 = float(pop3_)/float(t_pop)
rate4 = float(pop4_)/float(t_pop)
rate5 = float(pop5_)/float(t_pop)
rate6 = float(pop6_)/float(t_pop) 

print("CASE1: ", pop1_, "rate: ", rate1)
print('\n')
print("CASE2: ", pop2_, "rate: ", rate2)
print('\n')
print("CASE3: ", pop3_, "rate: ", rate3)
print('\n')
print("CASE4: ", pop4_, "rate: ", rate4)
print('\n')
print("CASE5: ", pop5_, "rate: ", rate5)
print('\n')
print("CASE6: ", pop6_, "rate: ", rate6)
print('\n')
print("total: ", t_pop)

print([pop1_/t_pop, pop2_/t_pop, pop3_/t_pop, pop4_/t_pop, pop5_/t_pop, pop6_/t_pop])

