#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=b1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6


  pwd=`pwd`

  for i in `seq 1 1 10`
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