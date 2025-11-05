#!/bin/sh
#SBATCH --partition=3080ti
#SBATCH --job-name=DB_SCF
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --gpus-per-task=1

module load compiler mkl mpi
module load cuda/11.6
module load pwmat

 mpirun -np $SLURM_NPROCS PWmat | tee output




