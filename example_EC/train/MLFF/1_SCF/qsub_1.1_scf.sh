#!/bin/bash

  # ----- setting parameter start ----- 

  T_list=("400K" "700K" "1200K")

  # ----- setting parameter end   ----- 

  echo " ###### start copy qsub ###### "

  pwd=`pwd`

  for it in 0 1 2
  do
    T=${T_list[$it]}

    echo " ----- ${T} ----- "

    rm -rf ${T} && mkdir ${T}
    cp ./template/qsub_${T}.sh ./${T}
    cd ${T}
    # bash qsub_${T}.sh
    cd ${pwd}
    echo " @@@ Please go to ${T} directry to change "
    echo "     the parameters and variable 'sid' in qsub_${T}.sh "
    echo "     to obtain node balance."
  done

