#!/bin/bash


  for i in `seq 0 200`
  do
    jId=$((${i} + 974266))
    scancel ${jId}
  done
