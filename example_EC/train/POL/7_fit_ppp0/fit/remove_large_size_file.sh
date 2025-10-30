#!/bin/bash

  pwd=`pwd`


  for i in `seq 1 1 30`
  do
    wdir=atom_${i}
    echo ${wdir}

    cd ${wdir}/temp_folder
    rm -rf feature movement

    cd fitting
    rm -rf 0 1 trainData.txt.Ftype*

    cd ${pwd}
  done


  for i in `seq 1 1 20`
  do
    wdir=bond_${i}
    echo ${wdir}

    cd ${wdir}/temp_folder
    rm -rf feature movement

    cd fitting
    rm -rf 0 1 trainData.txt.Ftype*
    
    cd ${pwd}
  done