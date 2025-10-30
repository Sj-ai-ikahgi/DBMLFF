#!/bin/bash

  pwd=`pwd`

  modName=${pwd}/ppp_model_neigh_3

  rm -rf ${modName}
  mkdir ${modName}

  cp -rp ./results_plot ${modName}

  rm -rf ${modName}/ppp_neighbore
  mkdir ${modName}/ppp_neighbore

  cp -p ../parameters/MOVEMENT.type ${modName}/ppp_neighbore
  cp -p ../parameters/IN.NEIGHBORE_3/IN.NEIGHBORE_all   ${modName}/ppp_neighbore/IN.NEIGHBORE


  mdir=${modName}/ppp_model
  rm -rf ${mdir}
  mkdir ${mdir}

  for i in `seq 1 1 30`
  do
    wdir=atom_${i}
    echo ${wdir}

    cd ${wdir}/temp_folder/fitting

    rm -rf ${mdir}/${wdir}
    mkdir ${mdir}/${wdir}
    cp -rp input output RMSE ${mdir}/${wdir}
    cd fread_dfeat
    cp -p feat_PV.* feat.info fit_linearMM.input linear_fitB.ntype ${mdir}/${wdir}

    cd ${pwd}
  done


  for i in `seq 1 1 20`
  do
    wdir=bond_${i}
    echo ${wdir}

    cd ${wdir}/temp_folder/fitting

    rm -rf ${mdir}/${wdir}
    mkdir ${mdir}/${wdir}
    cp -rp input output RMSE ${mdir}/${wdir}
    cd fread_dfeat
    cp -p feat_PV.* feat.info fit_linearMM.input linear_fitB.ntype ${mdir}/${wdir}
    
    cd ${pwd}
  done