#!/bin/bash

  for i in `seq 1 1 5`
  do
    sbatch run_${i}.sh
  done





