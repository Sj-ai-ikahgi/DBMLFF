#!/bin/sh
#SBATCH --partition=cpu 
#SBATCH --job-name=LDT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

	wdir=`pwd`


	sysName=("20LiFSI-24DME-60TTE")

	infile_format=("DBMLFF")

	odir=("data_out/20LiFSI-24DME-60TTE")

	for i in 0
	do

    sysId=${sysName[$i]}_${infile_format[$i]} 
		todir=${wdir}/data_out
		
		python3 msd_compose.py ${sysId} ${todir} 
		sleep 1
	done

