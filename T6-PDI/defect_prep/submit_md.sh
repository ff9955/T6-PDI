#!/bin/bash

#SBATCH --job-name=CP2K-MD
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --account=e05-biosoft-blu
#SBATCH --partition=standard
#SBATCH --qos=standard

# Setup the batch environment
#module restore /etc/cray-pe.d/PrgEnv-gnu
#module swap cray-mpich/8.0.16 cray-mpich-ucx/8.0.16  
#module swap craype-network-ofi craype-network-ucx 
module load cp2k/cp2k-2023.2

export OMP_NUM_THREADS=1


# (submit every 5 starting from frame 5 until 100)

srun --hint=nomultithread --distribution=block:block cp2k.psmp -i run.inp > run.log
