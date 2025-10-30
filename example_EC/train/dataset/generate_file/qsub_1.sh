#!/bin/bash


  # ----- setting parameter start ----- 

  whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
  export PYTHONPATH=${whome}/dataset:${PYTHONPATH}
  molType="EC"

  # ----- setting parameter end   ----- 


  echo "--- ${molType} --- "
  python generate_nearest_neigh_and_bond_len.py  ${molType}
  python generate_bond_dis_matrix.py             ${molType}
  python generate_in.neighbore.py                ${molType}