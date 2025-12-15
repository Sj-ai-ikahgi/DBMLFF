
all: cuda_routine

# ========== Dependency Library Paths ==========
# Use environment variables with default fallback values
MPI_HOME ?= /usr/local/mpich
NCCL_HOME ?= /usr/local/nccl

# ========== Compilation Flags ==========
# Define include and library paths using environment variables
MPI_INC = -I $(MPI_HOME)/include
MPI_LIB = -L $(MPI_HOME)/lib
NCCL_INC = -I $(NCCL_HOME)/include
NCCL_LIB = -L $(NCCL_HOME)/lib

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