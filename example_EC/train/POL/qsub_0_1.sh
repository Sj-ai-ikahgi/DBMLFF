#!/bin/bash

  # ----- setting parameter start ----- 

  partition="1080ti"

  # ----- setting parameter end   -----

  
  pwd=`pwd`

  # --- start first relax structure
  cd ./0_relax_structure/relax_1
  bash remove_pwmat.sh
  cp run_template.sh run.sh
  sed -i "s!@partition!${partition}!g" run.sh
  jobId=$(sbatch --parsable run.sh)
  echo "@@@ 1 ${jobId} "
  cd ${pwd}

  # --- start second relax structure
  cd ./0_relax_structure/relax_2
  bash remove_pwmat.sh
  cp run_template.sh run.sh
  sed -i "s!@partition!${partition}!g" run.sh
  jobId=$(sbatch --parsable --dependency=afterok:${jobId}:+1  run.sh  )
  echo "@@@ 2 ${jobId} "
  cd ${pwd}

  # --- start get energy Etot0
  cd ./1_SCF
  bash remove_pwmat.sh
  cp run_template.sh run.sh
  sed -i "s!@partition!${partition}!g" run.sh
  jobId=$(sbatch --parsable --dependency=afterok:${jobId}:+1  run.sh  )
  echo "@@@ 3 ${jobId} "