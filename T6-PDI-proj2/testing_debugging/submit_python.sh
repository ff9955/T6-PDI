#!/bin/bash

#SBATCH --job-name=cp2k
#SBATCH --nodes=1
#SBATCH --tasks-per-node=50
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --account=e05-biosoft-blu
#SBATCH --partition=standard
#SBATCH --qos=short
##SBATCH --reservation=shortqos

module load cray-python

python xsh_analysis_duplicate.py
