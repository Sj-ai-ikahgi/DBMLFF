#!/bin/bash


  for i in `seq 0 200`
  do
    jId=$((${i} + 972286))
    scancel ${jId}
  done
