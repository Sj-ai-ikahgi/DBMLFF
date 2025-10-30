#!/bin/bash

  cp -p ../1_SCF/xatom0.config    ./
  cp -p ../1_SCF/Etot0            ./
  cp -p ../3_POL_DFT/polarization.out_all ./polarization.out
  ./convert2pxyz.x


