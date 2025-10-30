#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=n1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

module load intel/2020
module load mkl

  # ifort rho_fitV33.f90 -mkl -o rho_fitV33.x

  # fitting rho
  cp ./rho_fitV33.input_1   ./rho_fitV33.input
  ./rho_fitV33.x

  mv ./funcr_atom.fit.bin   ./funcr_atom.fit.bin.in
  cp ./rho_fitV33.input_2   ./rho_fitV33.input
  ./rho_fitV33.x


  # check the fitting
  confNum=50
  cp ./data/xatom${confNum}.config ./rho2xsf/atom.config
  cp ./OUT.FIT_rho                 ./rho2xsf
  cd rho2xsf
  ./convert_local.x  OUT.FIT_rho
  cd ..


  # smooth the rho
  cp ./funcr_atom.fit.bin  ./Smoothen_funcr/funcr_atom.fit.bin_origin
  cd ./Smoothen_funcr
  ./funcr_fit_Ecut.x  > funcr_fit_Ecut_results.txt
  ntype_cent=`cat funcr_fit_Ecut_results.txt | grep "ntype_cent" | awk '{print $2}'`
  python plot_rho.py ${ntype_cent}




