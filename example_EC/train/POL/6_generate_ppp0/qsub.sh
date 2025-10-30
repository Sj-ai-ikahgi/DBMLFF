#!/bin/bash


  pwd=`pwd`


  # ----- setting parameter start ----- 

  whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
  mpichHome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC/mpich-3.4.3"
  partition="1080ti"
  sId_list=("400K" "700K" "1200K" "400K_s" "600K_s" "1900K_s")
  confNum_list=(3000 1500 1500 5000 5000 5000)


  # ----- setting parameter end   ----- 


  codePath="${whome}/dataset/md_src_ppp0"
  dataPath="${whome}/MLFF/MLFF_dataset"
 
  # set model parameter
  rm -rf ./template/mol.1  &&  mkdir ./template/mol.1

  # copy DB
  cp -p ${whome}/DB/2_fit_DB/Smoothen_funcr/funcr_atom.fit.bin_new \
   ./template/mol.1/funcr_atom.fit.bin

  # copy MLFF
  cp -p ${whome}/MLFF/MLFF_dataset/parameters/mol.info ./template/mol.1
  cd ${whome}/MLFF/3_fit_mlff/fit_dir
  cp -rp input output ${whome}/POL/6_generate_ppp0/template/mol.1
  cd fread_dfeat
  cp -p bond.molecule dist_pair.in feat_PV.* feat.info \
     fit_linearMM.input linear_fitB.ntype \
     ${whome}/POL/6_generate_ppp0/template/mol.1

  cd ${pwd}

  # copy POL - Dij
  cp -p ../POL_dataset/polar_param.input_ppp0 ./template/mol.1
  cp -p ../5_fit_direct_pxyz/fit.BB ./template/mol.1


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
    sed -i "s!@dataPath!${dataPath}!g"   run_mpich.sh
    sed -i "s!@mpichHome!${mpichHome}!g" run_mpich.sh
    sed -i "s!@confNum!${confNum}!g"     run_mpich.sh
    sed -i "s!@sId!${sId}!g"             run_mpich.sh

    bash remove.sh
    sbatch run_mpich.sh
    cd ..
  done