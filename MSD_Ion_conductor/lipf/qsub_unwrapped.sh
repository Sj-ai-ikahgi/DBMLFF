#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=gro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


	wdir=`pwd`

	sysName=("80EC-10LiPF6-60DMC")

	infile_format=("DBMLFF")

	for i in 0
	do

		cp run_unwrapped.sh run_unwrapped_${i}.sh

		todir=${wdir}/data_out
		JN=UW_$i

		# ddir="/data/home/sqjiang/fitting_smooth_rho_ion/check/MD/39EC_33EMC_51DMC_6VC_11LiPF6"
		# cat ${ddir}/MOVEMENT > ./MOVEMENT

		spoint=1
		epoint=10000000

		tidir="./"

		sed -i "s!@JN!${JN}!g"                             run_unwrapped_${i}.sh
		sed -i "s!@sysName!${sysName[$i]}!g"               run_unwrapped_${i}.sh
		sed -i "s!@infile_format!${infile_format[$i]}!g"   run_unwrapped_${i}.sh
		sed -i "s!@idir!${tidir}!g"                     	 run_unwrapped_${i}.sh
		sed -i "s!@odir!${todir}!g"                        run_unwrapped_${i}.sh
		sed -i "s!@spoint!${spoint}!g"                     run_unwrapped_${i}.sh
		sed -i "s!@epoint!${epoint}!g"                     run_unwrapped_${i}.sh
		sed -i "s!@out!unwrapped_out.${i}!g"               run_unwrapped_${i}.sh

		# sbatch run_unwrapped_${i}.sh
		bash run_unwrapped_${i}.sh

		sleep 1

	done

wait


