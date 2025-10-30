#!/bin/sh
#SBATCH --partition=1080ti
#SBATCH --job-name=energy_shift
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --gpus-per-task=1

module load compiler mkl mpi
module load cuda/11.6
module load pwmat

  cp ../../MLFF_dataset/Conf_400K/DFT_1.config ./atom.config

  mpirun -np $SLURM_NPROCS PWmat | tee output

  sleep 3

  python get_energy_shift.py


wait



