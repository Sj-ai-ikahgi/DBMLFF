#!/bin/sh
#SBATCH --partition=3080ti
#SBATCH --job-name=200K
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
MPICH_HOME=/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EA/mpich-3.4.3
export LD_LIBRARY_PATH="${MPICH_HOME}/lib:${LD_LIBRARY_PATH}"
export CPATH="${MPICH_HOME}/include:${CPATH}"
export PATH="${MPICH_HOME}/bin:${PATH}"

#nccl
NCCL_HOME="/share/app/nvhpc22.2/Linux_x86_64/22.2/comm_libs/nccl"
export LD_LIBRARY_PATH="${NCCL_HOME}/lib:${LD_LIBRARY_PATH}"


codePath=/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EA/dataset/md_src_ppp0

rm -f OUT.* main_MD.x_*

for cid in `seq 1 200`
do
  echo ${cid}
  cp ../Conf/new_structures/structure_${cid}.pwmat ./atom.config
  mpirun -np ${SLURM_NPROCS} ${codePath}/main_MD.x >> main_MD.x_out 2>&1 >> main_MD.x_err
done

exit 0


