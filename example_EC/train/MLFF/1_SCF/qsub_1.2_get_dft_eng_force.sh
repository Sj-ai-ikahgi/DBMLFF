#!/bin/bash

  # ----- setting parameter start ----- 

  natoms=10
  total_conf_list=(3000 1500 1500)
  ssize_list=(200 200 200)
  T_list=("400K" "700K" "1200K")

  # ----- setting parameter end   ----- 

  pwd=`pwd`


  echo " ###### start check dft missing ###### "
  for it in 0 1 2
  do
    T=${T_list[$it]}
    total_conf=${total_conf_list[$it]}
    ssize=${ssize_list[$it]}

    echo " ----- ${T} ----- "

    cp ./template/check_missing.py ./${T}
    cd ${T}
    python check_missing.py ${natoms} ${total_conf} ${ssize}
    cd ${pwd}
  done


  echo " ###### start get DFT energy and force ###### "
  for it in 0 1 2
  do
    T=${T_list[$it]}
    total_conf=${total_conf_list[$it]}
    ssize=${ssize_list[$it]}

    echo " ----- ${T} ----- "

    cp ./template/get_dft_eng_force.py ./${T}
    cd ${T}
    python get_dft_eng_force.py ${natoms} ${total_conf} ${ssize} 

    cd ${pwd}

    echo " ### the result is stored in ./${T}/dft_eng_froce.dat ### "
  done