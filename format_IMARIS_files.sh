#!/bin/bash

#Format IMARIS files
#IMARIS outputs a folder with many files for each metrics per analyzed images
#Consolidated tables build with R

mkdir table

for i in `ls -d *_Statistics`
do
	cd $i
	for j in `ls *.csv`
	do
		tail -n +5 $j | cut -d "," -f 1 > `basename $j .csv`.txt
		paste *.txt | column -s '\t' -t > ../table/`basename $i _Statistics`.tsv
	done
	rm *.txt
	cd ..
done