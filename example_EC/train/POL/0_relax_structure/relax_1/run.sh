#!/bin/sh
#SBATCH --partition=1080ti
#SBATCH --job-name=relax_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --gpus-per-task=1

module load compiler mkl mpi
module load cuda/11.6
module load pwmat

  cp -p ../atom.config ./atom.config
  mpirun -np $SLURM_NPROCS PWmat | tee output


