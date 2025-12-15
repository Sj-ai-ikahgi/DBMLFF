all: main_MD.x


# ========== Dependency Library Paths ==========
# Use environment variables with default fallback values
MPI_HOME ?= /usr/local/mpich
CUDA_HOME ?= /usr/local/cuda
NCCL_HOME ?= /usr/local/nccl
MKL_HOME ?= /opt/intel/mkl

# ========== Compilation Flags ==========
MPI_INC = -I $(MPI_HOME)/include
MPI_LIB = -L $(MPI_HOME)/lib
NCCL_INC = -I $(NCCL_HOME)/include
NCCL_LIB = -L $(NCCL_HOME)/lib -lnccl
MKL_INC = -I $(MKL_HOME)/include
MKL_LIB = -L $(MKL_HOME)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
CUDA_TOOLKIT_ROOT_DIR = $(CUDA_HOME)


MPIFC = pgfortran

FFLAG = -O3\
        -Mallocatable=95 -Mmpi=mpich -Mcuda \
        $(MPI_INC)  $(MPI_LIB) \
	      $(NCCL_INC) $(NCCL_LIB) -lnccl \
	    	$(MKL_INC)  $(MKL_LIB) -lstdc++ \
				-Mmkl \
		     -lpthread  -lm \
		    -I ./iotk/src/ \

%.o:  %.f90
	$(MPIFC) $(FFLAG) -c $<
%.o:  %.f
	$(MPIFC) $(FFLAG) -c $<


OBJ = \
	mod_cuinterface.o \
	mod_profile.o \
	mod_mem_preallocated.o \
  mod_matom.o \
	mod_param_densityFF.o \
	mod_mpi.o \
	mod_data.o \
	mod_control.o \
	mod_md.o \
	mod_ppp.o \
	common_module_99.o \
	calc_lin_forceM.o \
	calc_lin_forceM_ppp.o \
	calc_VV_forceM.o \
	calc_NN_forceM.o \
	find_neighbore.o \
	find_neighbore_ppp.o \
	readf_xatom_new.o \
	md.o \
	ranfx.o \
	rantest.o \
	numeric_container.o \
	get_ALI.o \
	gaussj.o \
	gen_2b_feature.o \
	gen_3b_feature.o \
	gen_2b_feature_ppp.o \
	gen_3b_feature_ppp.o \
	find_feature_2b_type3.o \
	find_feature_3b_type3.o \
	two_atom_distance.o \
	ML_FF_EF.o \
	ML_FF_PPP.o \
	nonbond_force.o \
	calc_polar_force.o \
	calc_dvxc.o \
	calc_coul.o \
	calc_coul_xc_kin_V.o \
	cfft.o \
	cfftd.o \
	combine_FF.o \
	UxcCA.o \
	main_MD.o \
	read_key_words.o \
	scan_key_words.o \
	transform_to_upper.o \
	mod_setting.o \
	mod_calc_pol_dir.o


mod_cuinterface.mod: mod_cuinterface.f90
	$(MPIFC) $(FFLAG) -c $^
mod_profile.mod: mod_profile.f90
	$(MPIFC) $(FFLAG) -c $^
mod_mem_preallocated.mod: mod_mem_preallocated.f90
	$(MPIFC) $(FFLAG) -c $^
mod_mpi.mod: mod_mpi.f90
	$(MPIFC) $(FFLAG) -c $^	
calc_ftype1.mod: gen_2b_feature.f90
	$(MPIFC) $(FFLAG) -c $^
calc_ftype2.mod: gen_3b_feature.f90
	$(MPIFC) $(FFLAG) -c $^
calc_ftype1_ppp.mod: gen_2b_feature_ppp.f90
	$(MPIFC) $(FFLAG) -c $^
calc_ftype2_ppp.mod: gen_3b_feature_ppp.f90
	$(MPIFC) $(FFLAG) -c $^
calc_lin.mod: calc_lin_forceM.f90
	$(MPIFC) $(FFLAG) -c $^
calc_lin_ppp.mod: calc_lin_forceM_ppp.f90
	$(MPIFC) $(FFLAG) -c $^
# calc_ftype2b_3b.mod: gen_feature_2b_3b.f90
# 	$(MPIFC) $(FFLAG) -c $^
# calc_lin_PC.mod: calc_lin_forceM_PC.f90
# 	$(MPIFC) $(FFLAG) -c $^	
calc_nn.mod: calc_NN_forceM.f90
	$(MPIFC) $(FFLAG) -c $^
calc_vv.mod: calc_VV_forceM.f90
	$(MPIFC) $(FFLAG) -c $^
common_module_99.mod: common_module_99.f90
	$(MPIFC) $(FFLAG) -c $^
mod_parameter.mod: mod_matom.f90
	$(MPIFC) $(FFLAG) -c $^
mod_control.mod: mod_control.f90
	$(MPIFC) $(FFLAG) -c $^
mod_data.mod: mod_data.f90
	$(MPIFC) $(FFLAG) -c $^
mod_md.mod: mod_md.f90 
	$(MPIFC) $(FFLAG) -c $^
mod_param_densityff.mod: mod_param_densityFF.f90
	$(MPIFC) $(FFLAG) -c $^
mod_ppp.mod: mod_ppp.f90
	$(MPIFC) $(FFLAG) -c $^
mod_calc_pol_dir.mod: mod_calc_pol_dir.f90
	$(MPIFC) $(FFLAG) -c $^
mod_setting.mod: mod_setting.f90
	$(MPIFC) $(FFLAG) -c $^
find_neighbore_ppp.mod: find_neighbore_ppp.f90
	$(MPIFC) $(FFLAG) -c $^

mods = mod_cuinterface.mod \
	mod_profile.mod \
	mod_mem_preallocated.mod \
	mod_mpi.mod \
	calc_ftype1.mod \
	calc_ftype2.mod \
	calc_lin.mod \
	calc_nn.mod \
	calc_vv.mod \
	common_module_99.mod \
	mod_parameter.mod \
	mod_control.mod \
	mod_data.mod \
	mod_md.mod \
	mod_param_densityff.mod \
	mod_calc_pol_dir.mod \
	find_neighbore_ppp.mod \
	calc_ftype1_ppp.mod \
	calc_ftype2_ppp.mod \
	mod_setting.mod \
	mod_ppp.mod \
	calc_lin_ppp.mod


main_MD.x: $(mods) $(OBJ) clean-x
	make -C ./cuda_subroutine/src
	$(MPIFC) -O3 -Mcuda=rdc -c mod_cuinterface.f90
	$(MPIFC)  $(FFLAG) \
	./cuda_subroutine/src/nonbond_cukernel.o \
	./cuda_subroutine/src/util.o \
	./cuda_subroutine/src/calc_coul_xc_kin_V_cukernel.o \
	./cuda_subroutine/src/calc_polar_force.o \
	./cuda_subroutine/src/calc_coul.o \
	./cuda_subroutine/src/calc_dvxc.o \
	./cuda_subroutine/src/init.o \
	$(OBJ) \
	./iotk/src/libiotk.a \
	-L ${CUDA_TOOLKIT_ROOT_DIR}/lib64 -lcudart -lcufft \
	-o $@

clean-x:
	rm -f main_MD.x

clean:
	rm -f *.o *.mod *.x

clean-all:
	rm -f *.o ./cuda_subroutine/src/*.o *.mod *.x