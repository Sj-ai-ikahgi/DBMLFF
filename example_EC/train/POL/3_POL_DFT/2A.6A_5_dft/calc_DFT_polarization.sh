#!/bin/bash
#SBATCH --partition=1080ti
#SBATCH --job-name=2A.6A_5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --gpus-per-task=1

module load compiler mkl mpi
module load cuda/11.6
module load pwmat

 for i in `seq 1 1 199`
 do
    echo ${i}
    ./gen_config2.x $i
    mpirun -np $SLURM_NPROCS PWmat | tee output
    grep "E_tot(eV)    =" REPORT | cat >> polarization.out
    sleep 3
 done


