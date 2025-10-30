#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=mlff_fitting
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12



  # ----- setting parameter start ----- 

  whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"

  natoms=10
  # atomic type is based on MolAtomSysList
  atomTypeNum=5
  # O_TypeNum=2
  # C_TypeNum=2
  # H_TypeNum=1
  maxNeighborNum=10

  # ----- setting parameter end   ----- 

  mlffCodePath="${whome}/dataset/mlff_src"

  bash cp_file2template.sh

  rm -rf fit_dir
  cp -rp template fit_dir

  for T in "400K" "700K" "1200K" "400K_s" "600K_s" "1900K_s"
  do
    cp -rp ./data_set/${T} ./fit_dir
    cp -rp ./fit_dir/MOVEMENT.type ./fit_dir/${T}
  done

  cd fit_dir

  echo " @@@ generate linear fit parameter "
  python generate_fit_linearMM.input.py  ${atomTypeNum} ${maxNeighborNum} ${natoms}

  cp ${mlffCodePath}/parameters_template.py ${mlffCodePath}/parameters.py

  sed -i "s!@natoms!${natoms}!g"                  ${mlffCodePath}/parameters.py 
  sed -i "s!@atomTypeNum!${atomTypeNum}!g"        ${mlffCodePath}/parameters.py 
  # sed -i "s!@O_TypeNum!${O_TypeNum}!g"            ${mlffCodePath}/parameters.py 
  # sed -i "s!@C_TypeNum!${C_TypeNum}!g"            ${mlffCodePath}/parameters.py 
  # sed -i "s!@H_TypeNum!${H_TypeNum}!g"            ${mlffCodePath}/parameters.py 
  sed -i "s!@maxNeighborNum!${maxNeighborNum}!g"  ${mlffCodePath}/parameters.py 

  echo " @@@ generate features and calculate PCA "
  python ${mlffCodePath}/main.py 

  echo " @@@ linear fitting "
  mv bond.molecule ./fread_dfeat/bond.molecule
  mv fit_linearMM.input ./fread_dfeat
  cd ./fread_dfeat
  make linear_fitB.ntype

  echo " @@@ using model on the traning set "
  cp dist_pair.out dist_pair.in
  make energyL.pred.1

  cd ../
  python plot_dataset.py