#!/bin/bash


  pwd=`pwd`

  echo "----- polarization of single module -----"
  cd ./6_generate_ppp0
  bash qsub.sh
echo "-----Post-job execution: bash ./6_generate_ppp0/set_job_array.sh-----"  
  cd ${pwd}

  echo "----- start fitting atoms and bonds-----"
  cd ./7_fit_ppp0/parameters
  bash cp.sh
  cd ./7_fit_ppp0/parameters/IN.NEIGHBORE_3
  python generate_in.neighbore_for_all_atom.py
  python generate_in.neighbore_for_each_atom.py
  cd ./7_fit_ppp0/fit
  bash qsub.sh
  cd ${pwd}


  