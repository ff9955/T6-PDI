import multiprocessing as mp
import os, shutil, sys
import glob
#if not os.path.isdir("foo"):   os.mkdir("foo")
#if not os.path.isdir("bar"):   os.mkdir("bar")
#with open("foo/script.sh", "w") as f:
#    f.write("echo 'bob'")
#with open("bar/script.sh", "w") as f:
#    f.write("echo 'bob'")

nworkers = int(sys.argv[1])

ROOT_DIR = os.getcwd()
def run_scripts(dir_):
    print(dir_)
    os.chdir(dir_)
    #os.chmod("script.sh", int("755", base=8))
    os.system("pwd")
    os.system("/scratch/sgiannini/flavoured-cptk/cp2k/exe/cp2k_local/cp2k.sopt run.inp > run.log")
    os.chdir(ROOT_DIR)
    print "\n"

directories = glob.glob("run-fssh-*") #list of all directories
p = mp.Pool(nworkers)
p.map(run_scripts, directories)

#shutil.rmtree("foo")
#shutil.rmtree("bar")
