#!/bin/sh
#SBATCH --partition=1080ti
#SBATCH --job-name=egg_box
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --gpus-per-task=1


UDA_HOME="/share/app/cuda/cuda-11.3"
export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"
export CPATH="${CUDA_HOME}/include:${CPATH}"
export PATH="${CUDA_HOME}/bin:${PATH}"

#nvhpc
export NVARCH="$(uname -s)_$(uname -m)"
export NVCOMPILERS="/share/app/nvhpc22.2"
export MANPATH="${NVCOMPILERS}/${NVARCH}/22.2/compilers/man:${MANPATH}"
export PATH="${NVCOMPILERS}/${NVARCH}/22.2/compilers/bin:${PATH}"

#mkl
. /share/app/intel2020u4/parallel_studio_xe_2020/psxevars.sh

#mpich
MPICH_HOME=/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/mpich-3.4.3
export LD_LIBRARY_PATH="${MPICH_HOME}/lib:${LD_LIBRARY_PATH}"
export CPATH="${MPICH_HOME}/include:${CPATH}"
export PATH="${MPICH_HOME}/bin:${PATH}"

#nccl
NCCL_HOME="/share/app/nvhpc22.2/Linux_x86_64/22.2/comm_libs/nccl"
export LD_LIBRARY_PATH="${NCCL_HOME}/lib:${LD_LIBRARY_PATH}"

module load python/3.8.3


  codePath=/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/dataset/md_src

  for n123 in  70 80 90 100 120 140 
  do
    cp MD.input_template MD.input
    sed -i "s!@nx!${n123}!g" MD.input
    sed -i "s!@ny!${n123}!g" MD.input
    sed -i "s!@nz!${n123}!g" MD.input
    for c in `seq 0 19`
    do
      echo "run ${n123} ${c}"
      cp ./shift_config/atom_${c}.config atom.config
      mpirun -np ${SLURM_NPROCS} ${codePath}/main_MD.x >> main_MD.x_out 2>&1 >> main_MD.x_err
    done
    odir=${n123}
    rm -rf ${odir} && mkdir ${odir}
    mv OUT.* ${odir}
  done

  sleep 3
  python plot_egg_box_energy.py
  python plot_egg_box_force.py

exit 0