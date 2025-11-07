#!/bin/bash

#SBATCH --job-name=cp2k
#SBATCH --nodes=1
#SBATCH --tasks-per-node=50
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --account=e864
#SBATCH --partition=standard
#SBATCH --qos=short
##SBATCH --reservation=shortqos

source /work/e864/e864/fivanovic_extra/numba_venv/bin/activate

export OMP_NUM_THREADS=1

for i in $(seq $1 $2);do

	        python ../analysis_scripts/xsh_hamiltonian_Erank.py ../physical_system_optical $i &
	done
	wait
