#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=ppp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

  pwd=`pwd`

  for i in 11
  do
    wdir=atom_${i}
    echo ${wdir}

    rm -rf ${wdir}
    cp -rp template ${wdir}
    cd ${wdir}
    python workflow_ppp_model_fitting_ht.py
    python plot_dataset.py atom_${i}
    cd ${pwd}
  done


  for i in 1
  do
    wdir=bond_${i}
    echo ${wdir}

    rm -rf ${wdir}
    cp -rp template ${wdir}
    cd ${wdir}
    python workflow_ppp_model_fitting_ht.py
    python plot_dataset.py bond_${i}
    cd ${pwd}
  done



