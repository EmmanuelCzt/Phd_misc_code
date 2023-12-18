#!/bin/bash

mkdir gffcompare

gtf_list=$1
ref_annot=$2

for i in `cat ${gtf_list}`
do
	gffcompare -R -r ${ref_annot} -M -N -T -o gffcompare/$i FEELnc\_$i/$i\_filtered.gtf.lncRNA.gtf
done 

# Using -R because we are interesting in discovery power so we cannot be as precise as we would like 
#Using STARalign conda env