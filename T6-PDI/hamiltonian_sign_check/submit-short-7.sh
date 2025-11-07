#!/bin/bash

#SBATCH --job-name=CP2K-POD
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --account=e05-biosoft-blu
#SBATCH --partition=standard
#SBATCH --qos=short

# Setup the batch environment
#module restore /etc/cray-pe.d/PrgEnv-gnu
#module swap cray-mpich/8.0.16 cray-mpich-ucx/8.0.16  
#module swap craype-network-ofi craype-network-ucx 
module load cp2k/cp2k-2023.2

export OMP_NUM_THREADS=1
#export input=DN4T_crys_DIMER_1_coff12.xyz_DZVP-GTH_PBE-GTH

srun --hint=nomultithread --distribution=block:block cp2k.psmp CRYS_PAIRS_interface-dimer-7.txt_DZVP-GTH_PBE-GTH.inp > CRYS_PAIRS_interface-dimer-7.txt_DZVP-GTH_PBE-GTH.out
