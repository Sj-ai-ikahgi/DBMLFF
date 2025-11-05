#!/bin/sh
#SBATCH --partition=cpu
#SBATCH --job-name=LDT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

	wdir=`pwd`


	sysName=("80EC-10LiPF6-60DMC")

	infile_format=("DBMLFF")

	odir=("data_out/80EC-10LiPF6-60DMC")

	for i in 0
	do

    sysId=${sysName[$i]}_${infile_format[$i]} 
		todir=${wdir}/data_out
		
		python3 msd_compose.py ${sysId} ${todir} 
		sleep 1
	done

