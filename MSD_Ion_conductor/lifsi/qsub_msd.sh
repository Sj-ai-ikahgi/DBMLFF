#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=msd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


	wdir=`pwd`


	sysName=("20LiFSI-24DME-60TTE")

	infile_format=("DBMLFF")

	odir=("data_out/20LiFSI-24DME-60TTE")

	for i in 0
	do

		cp run_msd.sh run_msd_${i}.sh

    sysId=${sysName[$i]}_${infile_format[$i]} 
		todir=${wdir}/data_out
		JN=MSD_$i

		sed -i "s!@JN!${JN}!g"           run_msd_${i}.sh
		sed -i "s!@sysId!${sysId}!g"     run_msd_${i}.sh
		sed -i "s!@odir!${todir}!g"      run_msd_${i}.sh
		sed -i "s!@out!msd_out.${i}!g"   run_msd_${i}.sh

		bash run_msd_${i}.sh
		sleep 1
	done

