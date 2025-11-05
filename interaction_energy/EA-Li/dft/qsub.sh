#!/bin/bash

  # ----- setting parameter start ----- 

  partition="1080ti"
  nodes=1
  confNum=200

  # ----- setting parameter end   ----- 


  ppDir="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/FSI/dataset/Pseudopotential"

  cd ./z

  cp -p ./DFT_SCF/etot.input_template   ./DFT_SCF/etot.input
  sed -i "s!@ppDir!${ppDir}!g"          ./DFT_SCF/etot.input

  cp  -p run.sh_template                run.sh
  sed -i "s!@partition!${partition}!g"  run.sh
  sed -i "s!@nodes!${nodes}!g"          run.sh
  sed -i "s!@confNum!${confNum}!g"      run.sh

  sbatch run.sh