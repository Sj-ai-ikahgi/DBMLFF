#!/bin/bash

  # ----- setting parameter start ----- 

  whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
  mpichHome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/mpich-3.4.3"
  partition="1080ti"
  sId_list=("400K" "700K" "1200K" "400K_s" "600K_s" "1900K_s")
 confNum_list=(3000 1500 1500 5000 5000 5000)

  # ----- setting parameter end   ----- 


  codePath="${whome}/dataset/md_src"


  # set the mol charge and mol.info
  rm -rf ./template/mol.1  &&  mkdir ./template/mol.1
  cp -p ${whome}/MLFF/MLFF_dataset/parameters/mol.info ./template/mol.1
  cp -p ${whome}/DB/2_fit_DB/Smoothen_funcr/funcr_atom.fit.bin_new \
   ./template/mol.1/funcr_atom.fit.bin


  for i in `seq 0 5`
  do
    sId=${sId_list[$i]}

    wdir=${sId}

    rm -rf ${wdir} && cp -rp template ${wdir}
    cd ${wdir}


    jobName=${sId}
    confNum=${confNum_list[${i}]}

    sed -i "s!@partition!${partition}!g" run_mpich.sh
    sed -i "s!@job!${jobName}!g"         run_mpich.sh
    sed -i "s!@codePath!${codePath}!g"   run_mpich.sh
    sed -i "s!@mpichHome!${mpichHome}!g" run_mpich.sh
    sed -i "s!@confNum!${confNum}!g"     run_mpich.sh
    sed -i "s!@sId!${sId}!g"             run_mpich.sh

    bash remove.sh
    sbatch run_mpich.sh
    cd ..
  done