#!/bin/bash


  # ----- setting parameter start ----- 

  whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EA"
  partition="1080ti"

  # ----- setting parameter end   -----

  ppDir=${whome}/dataset/Pseudopotential 

  for T in 400 700 1200
  do
    rm -rf ${T}K
    cp -rp template ${T}K
    cd ${T}K
    mv etot.input_template etot.input
    sed -i "s!@T!${T}!g"   etot.input
    sed -i "s!@ppDir!${ppDir}!g"   etot.input
    sed -i "s!@partition!${partition}!g" submit.sh
    sed -i "s!@jobId!${T}!g" submit.sh
    sbatch submit.sh
    cd ..
  done