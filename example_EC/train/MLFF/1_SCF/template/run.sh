#!/bin/sh
#SBATCH --partition=@partition
#SBATCH --job-name=@JN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --gpus-per-task=1

module load compiler mkl mpi
module load cuda/11.6
module load pwmat


whome=@whome
ppDir=${whome}/dataset/Pseudopotential

for mid in `seq @spoint 1 @epoint`
do
  rm -rf ${mid} && mkdir ${mid}
  cd ${mid}

  cp ${whome}/MLFF/1_SCF/template/DFT_SCF/* .
  mv ./etot.input_template ./etot.input
  sed -i "s!@ppDir!${ppDir}!g"     etot.input

  icn="DFT_${mid}.config"
  cp ../conf/${icn} ./atom.config

  mpirun -np $SLURM_NPROCS PWmat | tee output
  sleep 3
  cd ../
done

wait



