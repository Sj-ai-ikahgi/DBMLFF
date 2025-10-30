
CUDA_HOME=/share/app/cuda/cuda-11.3
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export CPATH=$CUDA_HOME/include:$CPATH
export PATH=$CUDA_HOME/bin:$PATH

export NVARCH="$(uname -s)_$(uname -m)"
export NVCOMPILERS=/share/app/nvhpc22.2
export MANPATH=$NVCOMPILERS/$NVARCH/22.2/compilers/man:$MANPATH
export PATH=$NVCOMPILERS/$NVARCH/22.2/compilers/bin:$PATH

. /share/app/intel2020u4/parallel_studio_xe_2020/psxevars.sh


