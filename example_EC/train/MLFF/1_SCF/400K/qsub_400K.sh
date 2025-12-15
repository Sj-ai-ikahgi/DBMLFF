#!/bin/bash

  # ----- setting parameter start ----- 

   whome="/data/home/xiongrui1/sabaisheng/Shenjie/fitting/EC"
   partition="1080ti"

  # ----- setting parameter end   ----- 

   pwd=`pwd`

   ddir="${whome}/MLFF/MLFF_dataset/Conf_400K"

   for sid in `seq 0 1 14` 
   do
      spoint=$((1+${sid}*200))
      epoint=$((${spoint} + 200 - 1))

      wdir=${spoint}_${epoint}

      echo "${sid} ${wdir}"

      rm -rf ${wdir} && mkdir ${wdir}

      cp -p ../template/run.sh ${wdir}/run_${sid}.sh
      cd ${wdir}

      rm -rf conf && mkdir conf
      for cid in `seq ${spoint} 1 ${epoint}`
      do
         cp -p ${ddir}/DFT_${cid}.config ./conf
      done

      sed -i "s!@partition!${partition}!g" run_${sid}.sh
      sed -i "s!@JN!${sid}!g"              run_${sid}.sh
      sed -i "s!@spoint!${spoint}!g"       run_${sid}.sh
      sed -i "s!@epoint!${epoint}!g"       run_${sid}.sh
      sed -i "s!@whome!${whome}!g"         run_${sid}.sh
      sbatch  run_${sid}.sh

      cd ${pwd}
   done