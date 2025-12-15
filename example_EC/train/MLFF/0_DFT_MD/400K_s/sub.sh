#!/bin/bash

#SBATCH --partition=xahdnormal
#SBATCH --job-name=pwmat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=dcu:4
#SBATCH --exclusive

ulimit -a
ulimit -s unlimited
ulimit -a

module purge

module load pwmat/20240223-magma-DCU2-intelmpi-2021

mpirun  -np $SLURM_NPROCS PWmat -host 10.13.5.14 50001 | tee output
