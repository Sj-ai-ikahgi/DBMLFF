#!/bin/bash

  # ----- setting parameter start ----- 

  whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
  partition="1080ti"
  nodes=1
  confNum=50

  # ----- setting parameter end   ----- 


  ppDir="${whome}/dataset/Pseudopotential"

  bash clean_1.sh
  cp -p ./DB_dataset/Conf/DFT_*.config  ./1_SCF/Conf
  cd ./1_SCF

  cp -p ./DFT_SCF/etot.input_template   ./DFT_SCF/etot.input
  sed -i "s!@ppDir!${ppDir}!g"          ./DFT_SCF/etot.input

  cp  -p run.sh_template                run.sh
  sed -i "s!@partition!${partition}!g"  run.sh
  sed -i "s!@nodes!${nodes}!g"          run.sh
  sed -i "s!@confNum!${confNum}!g"      run.sh

  sbatch run.sh