#!/bin/bash

  pwd=`pwd`

  rm -rf results_plot && mkdir results_plot

  for i in `seq 1 1 30`
  do
    wdir=atom_${i}
    echo ${wdir}

    if [ ! -d ${wdir} ]; then
      echo " ${wdir} not exist"
      cycle
    fi

    cp ./template/plot_dataset.py  ${wdir}
    cd ${wdir}
    python plot_dataset.py  atom_${i}
    mv ./eng_force_dataset.png  ../results_plot/eng_atom_${i}.png
    cd ${pwd}
  done


  for i in `seq 1 1 20`
  do
    wdir=bond_${i}
    echo ${wdir}

    if [ ! -d ${wdir} ]; then
      echo " ${wdir} not exist"
      continue
    fi

    cp ./template/plot_dataset.py  ${wdir}
    cd ${wdir}
    python plot_dataset.py  bond_${i}
    mv ./eng_force_dataset.png  ../results_plot/eng_bond_${i}.png
    cd ${pwd}
  done