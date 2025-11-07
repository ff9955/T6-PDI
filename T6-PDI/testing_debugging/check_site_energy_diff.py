import numpy as np
import glob

"""
This script analyze the site energy difference to check the goodness of the FF.
Just give number of diabatic states and the difference you want to calculate
"""

datadict = "./"
ndiab = 12

def _diff_calc(s1,s2, glob_list):
    state1 = glob_list[s1-1]
    state2 = glob_list[s2-1]
    
    mean_diff = np.mean(np.array(state2) - np.array(state1))
    return mean_diff*27.2114*1000

glob_list = [[] for i in range(ndiab)]
for dir_ in glob.glob(datadict + "run-fssh-*"):
    filename = dir_ + "/" + "run-pseudo-hamilt-1.xyz"
    with open(filename, "r") as f:
        for line in f.readlines():
            if "nadiab" in line:
                pass
            elif "time" in line:
                time_ = line.split()[2].replace(",", "")
                if time_ != "0":
                    break
            else:
                splitted = line.split() 
                if splitted[0] == splitted[1]:
                    glob_list[int(splitted[0])-1].append(float(splitted[2]))
print "GOT glob_list"


#for i in range(2,ndiab+1):
#    print "state", i, i-1, _diff_calc(i,i-1, glob_list)


print "state 2, 1 (1-5,1-4) labelec", _diff_calc(2,1, glob_list)
print "state 3, 1 (1-6,1-4) labelec", _diff_calc(3,1, glob_list)
print "state 3, 2 (1-6,1-5) zero", _diff_calc(3,2, glob_list)
print "state 4, 1 (2-4,1-4) labhole", _diff_calc(4,1, glob_list)
print "state 5, 1 (2-5,1-4) labelec+lhole", _diff_calc(5,1, glob_list)
print "state 4, 5 (2-5,2-4) labelec", _diff_calc(4,5, glob_list)
