#!/bin/bash


  # ----- setting parameter start ----- 

whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
  partition="1080ti"

  # ----- setting parameter end   ----- 


  ppDir=${whome}/dataset/Pseudopotential

  cd energy_shift

  bash remove_pwmat.sh

  cp etot.input_template  etot.input
  sed -i "s!@ppDir!${ppDir}!g"           etot.input
  cp run_template.sh                     run.sh
  sed -i "s!@partition!${partition}!g"   run.sh
  sbatch run.sh


  echo " ### The result will be in ./energy_shift/energy_shift.dat ### "