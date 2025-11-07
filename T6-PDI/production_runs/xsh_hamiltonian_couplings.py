import numpy as np
from file_parsers import *
import sys

target_directory = sys.argv[1]
trajectory_list = [int(sys.argv[2])]

def build(pseudoH_lines, number_diabats):

    hamiltonian = np.zeros((number_diabats, number_diabats))

    for element in pseudoH_lines:

       row = int(element[0]) - 1
       column = int(element[1]) - 1
       quantity = element[2]

       if row == column:
           quantity = quantity/2

       hamiltonian[row,column] = quantity

    hamiltonian = hamiltonian + hamiltonian.T

    return hamiltonian

number_excitons=20
number_diabats=420
CT_number = number_diabats - number_excitons

ee_couplings = []
hh_couplings = []
xt_ct_couplings = []
xt_ct_deltaE_list = []
xt_couplings = []

#couplings between second-closest neighbours
hh_secondary = []
ee_secondary = []
xt_secondary = []

xtct_secondary_h = []
xtct_secondary_xt = []

for trajectory in trajectory_list:

    with open(f'{target_directory}/run-fssh-{trajectory}/run-pseudo-hamilt-1.xyz') as H_file:

        H_file.readline()
        H_file.readline()
        parsed_H = []
        counter = 0

        for line in H_file:

            split_line = line.split()
            if len(split_line) == 3:
                parsed_H.append(split_line)

            else:

                if counter%10 != 0:

                    parsed_H = []
                    H_file.readline()
                    counter += 1
                    continue

                parsed_H = np.array(parsed_H, dtype=float)
                H = build(parsed_H, number_diabats)

                parsed_H = []
                H_file.readline()
                counter += 1

                #interating over every 2nd donor molecule
                for j in range(0,(CT_number),2*number_excitons):
                    donor_row = H[j]

                    first_coupling = j+20 #h-h coupling with nearest neighbour to left of 6T
                    hh_couplings.append(donor_row[first_coupling])

                    if j >= (CT_number-2*number_excitons):
                        break
                    else:
                        second_coupling = j+60 #h-h coupling with nearest neighbour to right of 6T
                        hh_couplings.append(donor_row[second_coupling])

                        #h-h coupling between 2nd-nearest neighbours
                        non_nn_coupling = j+40
                        hh_secondary.append(donor_row[non_nn_coupling])

                        #therefore iterating in increments of 40 as couplings of 6T molecules inbetween these indices already taken into account
                for j in range(0,number_excitons-1):
                    acceptor_row = H[j]

                    coupling = j+1
                    non_nn_coupling = j+2

                    #e-e coupling between 2nd nearest neighbours
                    ee_couplings.append(acceptor_row[coupling])

                    if j <= (number_excitons-3):  ee_secondary.append(acceptor_row[non_nn_coupling])

                xt_ct_coupling = H[360, 400]
                xt_ct_couplings.append(xt_ct_coupling)

                #secondary xt-ct coupling between 2nd hole from interface and interfacial exciton
                xtct_2nd_exciton = H[361,401]
                xtct_secondary_xt.append(xtct_2nd_exciton)

                #secondary xt-ct coupling between 2nd exciton from interface and interfacial hole
                xtct_2nd_hole = H[380,400]
                xtct_secondary_h.append(xtct_2nd_hole)

                xt_ct_deltaE = H[400,400] - H[360,360]
                xt_ct_deltaE_list.append(xt_ct_deltaE)

                xt_couplings.append(H[400,401])
                xt_couplings.append(H[406,407])

                #2nd-nearest neaighbour xt-xt couplings
                xt_secondary.append(H[400,402])
                xt_secondary.append(H[406,408])

filename = 'e5_2xCT'

data_file = open(f'{filename}_ee_couplings.txt', 'a')

for element in ee_couplings:
    data_file.write(str(element) + '\n')

data_file.close()

data_file = open(f'{filename}_hh_couplings.txt', 'a')

for element in hh_couplings:
    data_file.write(str(element) + '\n')

data_file.close()

data_file = open(f'{filename}_xtct_couplings.txt', 'a')
for element in xt_ct_couplings:
    data_file.write(str(element) + '\n')
data_file.close()

data_file = open(f'{filename}_xtct_deltaE.txt', 'a')
for element in xt_ct_deltaE_list:
    data_file.write(str(element) + '\n')
data_file.close()

data_file = open(f'{filename}_xt_couplings.txt', 'a')
for element in xt_couplings:
    data_file.write(str(element) + '\n')
data_file.close()

data_file = open(f'{filename}_xt_secondary_couplings.txt', 'a')
for element in xt_secondary:
    data_file.write(str(element) + '\n')
data_file.close()

data_file = open(f'{filename}_xtct_secondary_h_couplings.txt', 'a')
for element in xtct_secondary_h:
    data_file.write(str(element) + '\n')
data_file.close()

data_file = open(f'{filename}_xtct_secondary_xt_couplings.txt', 'a')
for element in xtct_secondary_xt:
    data_file.write(str(element) + '\n')
data_file.close()

data_file = open(f'{filename}_ee_secondary_couplings.txt', 'a')
for element in ee_secondary:
    data_file.write(str(element) + '\n')
data_file.close()

data_file = open(f'{filename}_hh_secondary_couplings.txt', 'a')
for element in hh_secondary:
    data_file.write(str(element) + '\n')
data_file.close()


