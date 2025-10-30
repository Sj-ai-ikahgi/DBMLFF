#!/bin/bash

  result_file=./polarization.out_all
  rm -f ${result_file}

  sys_list=@sys_list
  split_num_list=@split_num_list

  for si in `seq 0 1 5`
  do
    sys=${sys_list[$si]}
    split_num=${split_num_list[$si]}

    for p in `seq 1 1 ${split_num}`
    do
      wdir=${sys}_${p}
      echo ${wdir}

      # rm ${wdir}/OUT.WG
      # rm ${wdir}/OUT.GKK
      # rm ${wdir}/OUT.RHO
      # rm ${wdir}/OUT.VR
      # rm ${wdir}/OUT.VR_hion

      python combine_pxyz_eng.py  ${wdir} ${result_file}

      # python remove_DFT_not_calculate.py ${wdir}/polarization.out
      # cat temp.out >> ./polarization.out
    done
  done

  # rm temp.out
  # mv ./polarization.out ../convert2pxyz


