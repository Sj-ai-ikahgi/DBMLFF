#!/bin/sh
#SBATCH --partition=@partition
#SBATCH --job-name=SCF
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --gpus-per-task=1

module load compiler mkl mpi
module load cuda/11.6
module load pwmat

  cp -p ../0_relax_structure/relax_2/final.config ./xatom0.config
  python change_config_format.py

  mpirun -np $SLURM_NPROCS PWmat  | tee output
  grep "E_tot(eV)    =" REPORT | cat > Etot0

wait



