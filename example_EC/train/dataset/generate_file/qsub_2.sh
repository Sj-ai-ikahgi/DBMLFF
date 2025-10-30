#!/bin/bash


  # for molType in "EC" "PF6" "DEC" "DMC" "EMC" "PC" "VC" "EA"

  # ----- setting parameter start ----- 

  whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
  export PYTHONPATH=${whome}/dataset:${PYTHONPATH}
  molType="EC"

  # ----- setting parameter end   ----- 


  # --- generate parameter files
  echo "--- ${molType} --- "
  python generate_file.py              ${molType}
  python generate_atom_move_order.py   ${molType}


  # --- copy the parameter files to right directory
  cd ./results

  # copy parameter file for DB fitting
  ddir="${whome}/DB/DB_dataset/parameters"
  rm -rf ${ddir} && mkdir ${ddir}
  cp -p mol_atom_bond_info.dat_${molType}   ${ddir}/mol_atom_bond_info.dat
  cp -p mol.info_${molType}                 ${ddir}/mol.info

  # copy parameter files for MLFF fitting
  ddir="${whome}/MLFF/MLFF_dataset/parameters"
  rm -rf ${ddir} && mkdir ${ddir}
  cp -p move_atom_order.dat_${molType}   ${ddir}/move_atom_order.dat
  cp -p mol.info_${molType}              ${ddir}/mol.info
  cp -p find_neighbore.in                ${ddir}/find_neighbore.in
  cp -p IN.NEIGHBORE_${molType}          ${ddir}/IN.NEIGHBORE
  cp -p bond.molecule_${molType}_bt      ${ddir}/bond.molecule
  cp -p MOVEMENT.type_${molType}         ${ddir}/MOVEMENT.type