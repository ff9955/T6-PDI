#!/bin/bash

#SBATCH --job-name=cp2k
#SBATCH --nodes=1
#SBATCH --tasks-per-node=50
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --account=e864
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --reservation=shortqos

module load cray-python

export OMP_NUM_THREADS=1

for i in $(seq $1 $2)
do
   cd run-fssh-$i
       /work/e864/e864/fivanovic_extra/flavoured-cptk-X-SH-coulomb-barrier/cp2k/exe/archer2/cp2k.sopt run.inp > run.log &
   cd ../
done

wait

