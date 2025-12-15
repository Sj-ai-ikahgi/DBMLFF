#!/bin/bash
# ==============================================
# Environment Variables Configuration Template
# Please modify the following variables according to your system environment
# ==============================================

# ========== CUDA Configuration ==========
# Please set your CUDA installation path
export CUDA_HOME=/path/to/your/cuda  # e.g., /usr/local/cuda-11.3

# ========== MPI Configuration ==========
# Please set your MPI installation path
# Option 1: Use pre-compiled MPI in the project directory
# export MPI_HOME=$(pwd)/../../mpich-3.4.3
# Option 2: Use system MPI installation
export MPI_HOME=/path/to/your/mpi  # e.g., /usr/local/mpich

# ========== NVIDIA HPC SDK Configuration ==========
# Please set your NVIDIA HPC SDK installation path (for pgfortran compiler)
export NVCOMPILERS=/path/to/your/nvhpc  # e.g., /opt/nvidia/hpc_sdk

# ========== NCCL Configuration ==========
# Please set your NCCL installation path
export NCCL_HOME=/path/to/your/nccl  # e.g., /opt/nvidia/hpc_sdk/Linux_x86_64/22.2/comm_libs/nccl

# ========== Intel MKL Configuration ==========
# Please set your Intel MKL installation path
export MKL_HOME=/path/to/your/mkl  # e.g., /opt/intel/mkl
# pgfortran requires MKLROOT, not MKL_HOME
export MKLROOT=/path/to/your/mkl  # Same as MKL_HOME

# ========== Derived Environment Variables ==========
# The following variables are automatically built based on the above settings

# Set library paths
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${MPI_HOME}/lib:${NCCL_HOME}/lib:${MKL_HOME}/lib/intel64:${LD_LIBRARY_PATH}
# Note: MKL libraries are typically in lib/intel64, not just lib

# Set include paths
export CPATH=${CUDA_HOME}/include:${MPI_HOME}/include:${NCCL_HOME}/include:${MKL_HOME}/include:${CPATH}

# Set executable paths
export PATH=${CUDA_HOME}/bin:${MPI_HOME}/bin:${NVCOMPILERS}/Linux_x86_64/22.2/compilers/bin:${PATH}

# Set manual paths
export MANPATH=${NVCOMPILERS}/Linux_x86_64/22.2/compilers/man:${MANPATH}

# Load Intel compiler environment (optional, if available)
# Uncomment and modify the path if you have Intel compilers installed
# . /path/to/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64

echo "Environment variables configured successfully"
echo "CUDA_HOME: ${CUDA_HOME}"
echo "MPI_HOME: ${MPI_HOME}"
echo "MKL_HOME: ${MKL_HOME}"
echo "MKLROOT: ${MKLROOT}"
echo "NCCL_HOME: ${NCCL_HOME}"
echo "NVCOMPILERS: ${NVCOMPILERS}"