#!/bin/bash

  # ----- setting parameter start ----- 

  confNum=50
  z_valence=34.0
  # ----- setting parameter end   ----- 

  bash clean_2.sh

  cd ./2_fit_DB

  cp cp_rho_and_prepare_input.sh_template cp_rho_and_prepare_input.sh
  sed -i "s!@z_valence!${z_valence}!g"    cp_rho_and_prepare_input.sh

  # ----- copy rho and file
  cp -p ../DB_dataset/parameters/mol_atom_bond_info.dat ./
  bash cp_rho_and_prepare_input.sh  ${confNum}

  # ----- fitting rho
  cp  -p run.sh_template            run.sh
  sed -i "s!@confNum!${confNum}!g"  run.sh
  sbatch  run.sh