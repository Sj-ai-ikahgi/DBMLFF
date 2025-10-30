
all: cuda_routine

MPI_INC =  -I /data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/mpich-3.4.3/include
MPI_LIB =  -L /data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/mpich-3.4.3/lib
NCCL_INC = -I /share/app/nvhpc22.2/Linux_x86_64/22.2/comm_libs/nccl/include
NCCL_LIB = -L /share/app/nvhpc22.2/Linux_x86_64/22.2/comm_libs/nccl/lib

MPIFC = nvcc
#FFLAG = -D_GLIBCXX_USE_CXX11_ABI=0 -lstdc++ -lgfortran -dc \
#		    -arch=compute_80 -std=c++11 -lineinfo \
#        ${MPI_INC} ${MPI_LIB} ${NCCL_INC} ${NCCL_LIB}
FFLAG = -dc \
				-gencode=arch=compute_60,code=sm_60 \
				-gencode=arch=compute_61,code=sm_61 \
				-gencode=arch=compute_86,code=sm_86 \
        ${MPI_INC} ${MPI_LIB} \
				${NCCL_INC} ${NCCL_LIB} \
				-std=c++11 -lineinfo 

%.o:  %.cu
	$(MPIFC) $(FFLAG) -c $<

OBJ = \
	util.o \
	nonbond_cukernel.o \
	calc_polar_force.o \
	calc_coul_xc_kin_V_cukernel.o \
	calc_coul.o \
	calc_dvxc.o

cuda_routine: $(OBJ)
	$(MPIFC) $(FFLAG) init.cpp -o init.o

clean:
	rm -f ./*.o