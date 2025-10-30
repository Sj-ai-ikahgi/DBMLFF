#!/bin/bash


  # ----- setting parameter start -----   

  cp -p ../POL_dataset/polar_param.input ./polar_param.input

  # ----- setting parameter end   -----


  cp -p ../4_convert2pxyz/pxyz.outC ./pxyz.out
  ./fit_direct_pxyz.x
  python plot_test.fit.py
