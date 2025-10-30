#!/bin/bash


  # ----- setting parameter start ----- 

  partition="1080ti"

  cp -p ../POL_dataset/polar_param.input ./template_dft/polar_param.input

  # ----- setting parameter end   -----

  cp -p ../1_SCF/xatom0.config  ./template_dft
  cp -p ../1_SCF/Etot0          ./template_dft
  cp -p ../1_SCF/OUT.RHO        ./template_dft/OUT.RHO0

  pwd=`pwd`

  sys_list=@sys_list
  split_num_list=@split_num_list
  confNum_list=@confNum_list


  for si in `seq 0 1 5`
  do
    sys=${sys_list[$si]}
    split_num=${split_num_list[$si]}
    confNum=${confNum_list[${si}]}

    for p in `seq 1 1 ${split_num}`
    do
      wdir=${sys}_${p}_dft
      echo "qsub $wdir" 

      rm -rf ${wdir}
      cp -rp template_dft   ${wdir}
      cp ./split_point_result/point.${sys}.${p} ${wdir}/point.all

      cd ${wdir}

      echo ${p} ${confNum}

      sed -i "s!@partition!${partition}!g"    calc_DFT_polarization.sh
      sed -i "s!@jn!${sys}_${p}!g"            calc_DFT_polarization.sh
      sed -i "s!@confNum!${confNum}!g"        calc_DFT_polarization.sh

      sbatch calc_DFT_polarization.sh
      cd ${pwd}
    done
  done



