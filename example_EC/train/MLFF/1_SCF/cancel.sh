#!/bin/bash

  for i in `seq 0 100`
  do
    jId=$((${i}+965975))
    scancel ${jId}
  done

