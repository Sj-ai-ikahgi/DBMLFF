#!/bin/bash

  # ----- setting parameter start ----- 

  whome="/path/to/your/example_EC/train"   #Please use absolute paths
  mpichHome="/path/to/your/mpich"
  partition="3080ti"

  # ----- setting parameter end   ----- 

  codePath="${whome}/dataset/md_src_ppp0"
  dataPath="${whome}/MLFF/0_DFT_MD/400K"

  # set real space grid mesh
  nline=`grep "N123      =" ${dataPath}/REPORT`
  nx=`echo ${nline} | awk '{print $3}' `
  ny=`echo ${nline} | awk '{print $4}' `
  nz=`echo ${nline} | awk '{print $5}' `
  echo "Set nx=$nx ny=${ny} nz=${nz}"

  # ## set the mol charge and mol.info
  rm -rf ./template/mol.1  &&  mkdir ./template/mol.1
  cp -p ${whome}/MLFF/MLFF_dataset/parameters/mol.info                     ./template/mol.1
  cp -p ${whome}/DB/2_fit_DB/Smoothen_funcr/funcr_atom.fit.bin_new \
   ./template/mol.1/funcr_atom.fit.bin

  cp -rp ${whome}/MLFF/3_fit_mlff/fit_dir/input                            ./template/mol.1
  cp -rp ${whome}/MLFF/3_fit_mlff/fit_dir/output                           ./template/mol.1
  cp -rp ${whome}/MLFF/3_fit_mlff/fit_dir/fread_dfeat/bond.molecule        ./template/mol.1
  cp -rp ${whome}/MLFF/3_fit_mlff/fit_dir/fread_dfeat/dist_pair.in         ./template/mol.1
  cp -rp ${whome}/MLFF/3_fit_mlff/fit_dir/fread_dfeat/feat_PV.*            ./template/mol.1
  cp -rp ${whome}/MLFF/3_fit_mlff/fit_dir/fread_dfeat/feat.info            ./template/mol.1
  cp -rp ${whome}/MLFF/3_fit_mlff/fit_dir/fread_dfeat/fit_linearMM.input   ./template/mol.1
  cp -rp ${whome}/MLFF/3_fit_mlff/fit_dir/fread_dfeat/linear_fitB.ntype    ./template/mol.1
 
  cp -rp ${whome}/POL/5_fit_direct_pxyz/fit.BB ./template/mol.1

  cp -rp ${whome}/POL/7_fit_ppp0/fit/ppp_model_neigh_3/ppp_model    ./template/mol.1
  cp -rp ${whome}/POL/7_fit_ppp0/fit/ppp_model_neigh_3/ppp_neighbore ./template/mol.1
  cp -rp ${whome}/POL/7_fit_ppp0/parameters/polar_param.input_ppp0           ./template/mol.1


  wdir=md_test

  rm -rf ${wdir} && cp -rp template ${wdir}
  cd ${wdir}
  sed -i "s!@nx!${nx}!g"         MD.input
  sed -i "s!@ny!${ny}!g"         MD.input
  sed -i "s!@nz!${nz}!g"         MD.input

  jobName=md_test

  sed -i "s!@partition!${partition}!g" run_mpich.sh
  sed -i "s!@job!${jobName}!g"         run_mpich.sh
  sed -i "s!@mpichHome!${mpichHome}!g" run_mpich.sh
  sed -i "s!@codePath!${codePath}!g"   run_mpich.sh

  bash remove.sh
  #sbatch run_mpich.sh