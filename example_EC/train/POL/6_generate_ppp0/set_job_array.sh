#!/bin/bash

  ### ref:https://unix.stackexchange.com/questions/498517/how-to-find-filenames-sort-it-alphabetically-then-only-head-the-first-3-charact

  # ----- setting parameter start ----- 

  sId_list=("400K" "700K" "1200K" "400K_s" "600K_s" "1900K_s")

  # ----- setting parameter end   ----- 


  pwd=`pwd`

  result_file=${pwd}/job_array.txt
  rm -f ${result_file}

  for i in `seq 0 5`
  do
    T=${sId_list[$i]}
    cd ${T}/dataset

    dlist=`find -type d -printf '%f\n' | sort -n`
    fId=0
    for tdir in ${dlist}
    do
      if [ ${fId} -eq 0 ]; then
        fId=$((${fId} + 1))
      else
        echo ${pwd}/${T}/dataset/${tdir} >> ${result_file}
      fi
    done

    cd ${pwd}
  done
