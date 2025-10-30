#!/bin/bash
#SBATCH --partition=@partition
#SBATCH --job-name=@jobId
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --gpus-per-task=1

module load compiler mkl mpi
module load cuda/11.6
module load pwmat

mpirun -np $SLURM_NPROCS PWmat > pwmat_out 2> pwmat_err