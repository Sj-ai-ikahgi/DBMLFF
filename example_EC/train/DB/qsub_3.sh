#!/bin/bash

  # ----- setting parameter start ----- 

  whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
  mpichHome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/mpich-3.4.3"
  partition="1080ti"

  # ----- setting parameter end   ----- 

  export PYTHONPATH="${whome}/dataset":${PYTHONPATH}
  codePath="${whome}/dataset/md_src"

  bash clean_3.sh

  cd 3_check/egg_box
  cp -rp template ./sample
  cd ./sample

  cp -p ${whome}/DB/1_SCF/Conf/DFT_1.config ./atom.config
  python shift_mol.py

  sed -i "s!@partition!${partition}!g"  run_mpich.sh
  sed -i "s!@mpichHome!${mpichHome}!g"  run_mpich.sh
  sed -i "s!@codePath!${codePath}!g"    run_mpich.sh

  # set the mol charge and mol.info
  mkdir ./mol.1
  cp -p ${whome}/DB/DB_dataset/parameters/mol.info ./mol.1
  cp -p ${whome}/DB/2_fit_DB/Smoothen_funcr/funcr_atom.fit.bin_new \
   ./mol.1/funcr_atom.fit.bin
   
  sbatch run_mpich.sh