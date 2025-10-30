#!/bin/bash

  dhead="../../MLFF/MLFF_dataset"
  odir=./Conf

  rm -rf ${odir} && mkdir ${odir}


  RANDOM=2095
  cid=0

  spoint=1
  epoint=1200
  for i in `seq 1 20`
  do
    sid=`echo $(( RANDOM % (${epoint} - ${spoint} + 1 ) + ${spoint} ))`
    echo $i $sid
    cid=$((${cid} + 1))
    cp ${dhead}/Conf_400K/DFT_${sid}.config  ${odir}/DFT_${cid}.config
  done


  spoint=1
  epoint=500
  for i in `seq 1 15`
  do
    sid=`echo $(( RANDOM % (${epoint} - ${spoint} + 1 ) + ${spoint} ))`
    echo $i $sid
    cid=$((${cid} + 1))
    cp ${dhead}/Conf_700K/DFT_${sid}.config  ${odir}/DFT_${cid}.config
  done


  spoint=1
  epoint=500
  for i in `seq 1 15`
  do
    sid=`echo $(( RANDOM % (${epoint} - ${spoint} + 1 ) + ${spoint} ))`
    echo $i $sid
    cid=$((${cid} + 1))
    cp ${dhead}/Conf_1200K/DFT_${sid}.config  ${odir}/DFT_${cid}.config
  done