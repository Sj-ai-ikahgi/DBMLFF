#!/bin/sh
#SBATCH --partition=@partition
#SBATCH --job-name=@job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --gpus-per-task=1

# ========== CUDA Configuration ==========
# Please set your CUDA installation path
export CUDA_HOME=/path/to/your/cuda  # e.g., /usr/local/cuda-11.3
export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"
export CPATH="${CUDA_HOME}/include:${CPATH}"
export PATH="${CUDA_HOME}/bin:${PATH}"

# ========== NVIDIA HPC SDK Configuration ==========
# Please set your NVIDIA HPC SDK installation path
export NVARCH="$(uname -s)_$(uname -m)"
export NVCOMPILERS=/path/to/your/nvhpc  # e.g., /opt/nvidia/hpc_sdk
export MANPATH="${NVCOMPILERS}/${NVARCH}/22.2/compilers/man:${MANPATH}"
export PATH="${NVCOMPILERS}/${NVARCH}/22.2/compilers/bin:${PATH}"

# ========== Intel MKL Configuration ==========
# Please set your Intel MKL installation path
# Uncomment and modify the path if you have Intel compilers installed
 . /path/to/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64

# ========== MPI Configuration ==========
# Please set your MPI installation path
# This will be replaced by @mpichHome during job submission
MPICH_HOME=@mpichHome
export LD_LIBRARY_PATH="${MPICH_HOME}/lib:${LD_LIBRARY_PATH}"
export CPATH="${MPICH_HOME}/include:${CPATH}"
export PATH="${MPICH_HOME}/bin:${PATH}"

# ========== NCCL Configuration ==========
# Please set your NCCL installation path
export NCCL_HOME=/path/to/your/nccl  # e.g., /opt/nvidia/hpc_sdk/Linux_x86_64/22.2/comm_libs/nccl
export LD_LIBRARY_PATH="${NCCL_HOME}/lib:${LD_LIBRARY_PATH}"

# ========== Code Path Configuration ==========
# This will be replaced by @codePath during job submission
codePath=@codePath

# ========== Execution ==========
rm -f OUT.* main_MD.x_*
cp ../../train/MLFF/0_DFT_MD/template/atom.config ./atom.config
mpirun -np ${SLURM_NPROCS} ${codePath}/main_MD.x >> main_MD.x_out 2>&1 >> main_MD.x_err

exit 0