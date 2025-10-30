#!/bin/bash


  # ----- setting parameter start ----- 

  partition="cpu"

  cp -p ../POL_dataset/polar_param.input ./template_pxyz/polar_param.input

  # ----- setting parameter end   -----

  cp -p ../1_SCF/xatom0.config  ./template_pxyz
  cp -p ../1_SCF/Etot0          ./template_pxyz
  cp -p ../1_SCF/OUT.RHO        ./template_pxyz/OUT.RHO0

  pwd=`pwd`

  sys_list=(2A 4A 6A 2A.2A 2A.4A 2A.6A )
  split_num_list=(3 3 3 4 4 5 )
  confNum_list=(182 182 182 211 233 199 )


  for si in `seq 0 1 5`
  do
    sys=${sys_list[$si]}
    split_num=${split_num_list[$si]}
    confNum=${confNum_list[${si}]}

    for p in `seq 1 1 ${split_num}`
    do
      wdir=${sys}_${p}_pxyz
      echo "qsub $wdir" 

      rm -rf ${wdir}
      cp -rp template_pxyz   ${wdir}
      cp ./split_point_result/point.${sys}.${p} ${wdir}/point.all

      cd ${wdir}

      echo ${p} ${confNum}

      sed -i "s!@partition!${partition}!g"    calc_pxyz.sh
      sed -i "s!@jn!${sys}_${p}!g"            calc_pxyz.sh
      sed -i "s!@confNum!${confNum}!g"        calc_pxyz.sh

      sbatch calc_pxyz.sh
      cd ${pwd}
    done
  done



