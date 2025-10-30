#!/bin/bash

  export PYTHONPATH="../../../dataset":${PYTHONPATH}

  # ifort gen_points.f90 get_ALI.f gaussj.f -o gen_points.x
  # ifort choose_two.f90 ran1.f -o choose_two.x
  
  cp -p ../1_SCF/xatom0.config ./xatom0.config
  rm ./point.* ./charge_point_position/charge_position_*


  # generate one point file
  echo -e "2.0" | ./gen_points.x
  mv point.temp point.2A
  python remove_overlap_charge_position.py point.2A
  mv point.2A_remove point.2A
  p2ANum=`wc -l point.2A | awk '{print $1}'`
  p2ANum=$((${p2ANum}/2))

  echo -e "4.0" | ./gen_points.x
  mv point.temp point.4A
  python remove_overlap_charge_position.py point.4A
  mv point.4A_remove point.4A
  p4ANum=`wc -l point.4A | awk '{print $1}'`
  p4ANum=$((${p4ANum}/2))

  echo -e "6.0" | ./gen_points.x
  mv point.temp point.6A
  python remove_overlap_charge_position.py point.6A
  mv point.6A_remove point.6A
  p6ANum=`wc -l point.6A | awk '{print $1}'`
  p6ANum=$((${p6ANum}/2))


  # generate two point file
  echo -e "point.2A point.2A\n 125 ${p2ANum} ${p2ANum} " | ./choose_two.x
  mv point.all.2pt point.2A.2A

  echo -e "point.2A point.4A\n 135 ${p2ANum} ${p4ANum} " | ./choose_two.x
  mv point.all.2pt point.2A.4A

  echo -e "point.2A point.6A\n 145 ${p2ANum} ${p6ANum} " | ./choose_two.x
  mv point.all.2pt point.2A.6A


  # check charge point position
  cd ./charge_point_position
  python charge_point_position.py
  cd ..


