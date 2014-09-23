#!/bin/bash

# CDO algorithm
#---------------------------------------------------------------------------

dir_in=/home/nicholat/project/access/slabosc/
dir_out=/home/nicholat/project/access/slabosc/
file_list=$(<list.txt)

for i in $file_list
	do
	echo preparing $i

	newfile=$dir_in$i
	outfile=$dir_out$i
	echo newfile is $newfile 
	echo outfile is $outfile

# Get annual mean 
	cdo yearmean $outfile $outfile.ym.nc
# Detrend annual data
	cdo detrend $outfile.ym.nc $outfile.ym.dt.nc
# Remove seasonal cycle - rsc
#	cdo ymonsub $outfile.nt.nc -ymonavg $outfile.nt.nc $outfile.rsc.nc 
#	cdo detrend $outfile.rsc.nc $outfile.rsc.dt.nc



done


