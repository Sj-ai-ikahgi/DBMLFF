#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=4A_2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

 for i in `seq 1 1 182`
 do
    echo ${i}
    ./gen_config2.x $i
 done


