#!/bin/bash


  pwd=`pwd`

  echo "----- start generate pxyz.outC -----"
  cd ./4_convert2pxyz
  bash run.sh
  cd ${pwd}

  echo "----- start fitting -----"
  cd ./5_fit_direct_pxyz
  bash run.sh
  cd ${pwd}


  