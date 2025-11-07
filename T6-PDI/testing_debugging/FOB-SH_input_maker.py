

import os
import numpy as np 
import random
import glob
import re


def _clean_velocities( fileinname, fileoutname):
    filein = open(fileinname)
    fileout = open(fileoutname, 'w')
    for line in filein.readlines()[2:]:
        fileout.write( '\t'.join(map(str, line.split()[1:])) + '\n')
    filein.close()
    fileout.close()
    
def _cut_in_chunks(natoms_tot,input_lines, filenamein, filenameout):
    
    #crete additional directory 
    os.system("mkdir %s" %filenameout)
    
    lines_per_file = (natoms_tot) + input_lines
    smallfile = None
    count_frames = 0
    with open(filenamein) as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0:
                count_frames += 1
                if smallfile:
                    smallfile.close()
                small_filename = filenameout + "/"+ 'frame_%s.xyz' % (count_frames)
                smallfile = open(small_filename, "w")
            smallfile.write(line)
        if smallfile:
            smallfile.close()
    return count_frames
    

def _copy_positions_velocities(velocitities_position, repetition, pv_dir, datadir, pos=True):

    #starting from the end taking pos and vel
    list_indeces = []
    filename = pv_dir + '/frame_*'
    for elem in glob.glob(filename):
        indx = int(re.findall(r'\d+',elem)[-1])
        list_indeces.append(indx)
   
    print list_indeces 
    count = -1 
    for item in sorted(list_indeces)[-velocitities_position:]:# start from end
        pv =  pv_dir + '/' + 'frame_%s.xyz' % (item) #change this if you want last part of data
        for repeat in range(repetition):
            count +=1
            fob_sh_dir = datadir + '/' + "run-fssh-%s" % (count)
            
            print "COUNT", count
            print fob_sh_dir
            print pv
            os.system("cp %s %s" %(pv, fob_sh_dir))
            
            if pos is True:
                os.system("mv %s/frame* %s/pos-init.xyz" %(fob_sh_dir, fob_sh_dir))
            else:
                os.system("mv %s/frame* %s/vel-init.xyz" %(fob_sh_dir, fob_sh_dir))
                filenamein = "%s/vel-init.xyz" % fob_sh_dir
                filenameout = "%s/VELOC.init" % fob_sh_dir
                #create VELOC.init
                _clean_velocities(filenamein, filenameout)
                #os.system("rm %s/vel-init.xyz" % fob_sh_dir)
                
def _modify_inp(traj_to_make, datadir):
    for item in range(traj_to_make):
        fobsh_input = "%s/run-fssh-%s/run.inp" %( datadir, item)
        rand_num =  random.randint(1,1000000000)
        
        print fobsh_input
        print rand_num
        
        # Read in the file
        with open(fobsh_input, 'r') as file :
            filedata = file.read()
        
        # Replace the target string
        filedata = filedata.replace('SEED                                  2000', 'SEED                        %s' %rand_num)
        
        # Write the file out again
        with open(fobsh_input, 'w') as file:
            file.write(filedata)




############# MAIN ########################


#datadir = "/mnt/c/Users/samue/Desktop/Research_work/PENTACENE_THIN_FILMS/NEW_ADJUSTED_LAYERS/PREPARE_FOBSH/"
datadir = os.getcwd()

position_total  = datadir + "/run-pos-1.xyz"
velocity_total = datadir + "/run-vel-1.xyz"
pos_dir = datadir + "/pos"
vel_dir = datadir + "/vel"

#input cut chunks
#natoms_per_mol = 58152
#nmols = 1 
natoms_tot = 22800 
input_lines = 2

#run-fssh-* directories to create
velocitities_position = 50
repetition = 1
traj_to_make = velocitities_position*repetition

#create position frames
print "#######  CREATE POSITION FRAMES #######"
if not os.path.isdir(pos_dir): 
   count_frames = _cut_in_chunks(natoms_tot, input_lines, position_total, pos_dir)
else:
   print "POSITION DIRECTORY ALREADY PRESENT: USE THAT!"

#create velocity frames 
print "#######  CREATE VELOCITY FRAMES #######"
if not os.path.isdir(vel_dir):
   count_frames = _cut_in_chunks(natoms_tot, input_lines, velocity_total, vel_dir)
else:
   print "VELOCITY DIRECTORY ALREADY PRESENT: USE THAT!"

#copy template directory
for item in range(traj_to_make):
    os.system("cp -r %s/TO_COPY %s/run-fssh-%s" %(datadir, datadir, item))
    
# transfer position files
print "#######  TRANSFER POSITION FRAMES #######"
_copy_positions_velocities(velocitities_position, repetition, pos_dir, datadir, True)

# transfer velocity files 
print "#######  TRANSFER VELOCITY FRAMES #######"
_copy_positions_velocities(velocitities_position, repetition, vel_dir, datadir, False)

#modify input 
print "#######  MODIFY INPUT #######"
_modify_inp(traj_to_make, datadir)

print "DONE!"


