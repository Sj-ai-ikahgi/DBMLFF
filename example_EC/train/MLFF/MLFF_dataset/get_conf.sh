#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=get_conf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

  module load python/3.8.3


  sId_list=("400K" "700K" "1200K" "400K_s" "600K_s" "1900K_s")
  infile_format="DFT"
  molSize=10
  molNum_list=(15 15 15 1 1 1)
  confDelt_list=(20 40 40 1 1 1)

for i in `seq 3 5`
  do
    sId=${sId_list[$i]}
    idir="../0_DFT_MD/${sId}"
    odir="./Conf_${sId}"
    molNum=${molNum_list[$i]}
    confDelt=${confDelt_list[$i]}

    echo " --- start ${sId} --- "
    python get_conf.py ${sId} ${idir} ${odir} ${molNum} ${molSize} ${confDelt}
  done

  echo " ### check the result in ./temp ### "

# 运行统计脚本
echo "Running count_removed_molecules.py ..."
python count_removed_molecules.py

echo "统计完成！"