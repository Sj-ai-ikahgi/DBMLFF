#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=gro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

	wdir=`pwd`

	sysName=("20LiFSI-24DME-60TTE")

	infile_format=("DBMLFF")

	for i in 0
	do

		cp run_unwrapped.sh run_unwrapped_${i}.sh

		todir=${wdir}/data_out
		JN=UW_$i


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


