#!/bin/bash


  pwd=`pwd`

  echo "----- start generate pxyz.outC -----"
  cd ./6_generate_ppp0
  bash qsub.sh
  cd ${pwd}

  echo "----- start fitting -----"
  cd ./5_fit_direct_pxyz
  bash run.sh
  cd ${pwd}


  