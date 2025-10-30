#!/bin/bash
#SBATCH --partition=@partition
#SBATCH --job-name=@jn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

 for i in `seq 1 1 @confNum`
 do
    echo ${i}
    ./gen_config2.x $i
 done


