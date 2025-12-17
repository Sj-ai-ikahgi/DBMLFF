#!/bin/bash
#PBS -N pwmat
#PBS -l nodes=1:ppn=4
#PBS -q batch
#PBS -l walltime=20000:00:00

NPROCS=`wc -l < $PBS_NODEFILE`
cd $PBS_O_WORKDIR
ulimit -s unlimited
#export OMP_NUM_THREADS=24
#export OMP_STACKSIZE=512M
mpirun -np ${NPROCS} PWmat |tee pwmat.log

