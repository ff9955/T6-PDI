#!/bin/bash

#SBATCH --job-name=cp2k
#SBATCH --nodes=1
#SBATCH --tasks-per-node=50
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=e864
#SBATCH --partition=highmem
#SBATCH --qos=highmem

./store_checkpoint.sh $1 $2 &

export OMP_NUM_THREADS=1

for i in $(seq $1 $2)
do
   cd run-fssh-$i
       /work/e864/e864/fivanovic_extra/flavoured-cptk-X-SH-coulomb-barrier/cp2k/exe/archer2/cp2k.sopt run.inp > run.log &
   cd ../
done

wait

