#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=@JN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load intel/2020 
module load python/3.8.3

	python3 msd.py @sysId @odir | tee @out

