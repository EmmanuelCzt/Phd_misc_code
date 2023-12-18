#!/bin/bash


sample_list=$1
ref_annot=$2
folder=$3

mkdir ${folder}

for i in `cat ${sample_list}`
do
	mkdir ${folder}/$i
	/home/emmanuel/software/gffcompare-0.11.5.Linux_x86_64/gffcompare -o ${folder}/$i/$i\_R -r ${ref_annot} -R -D $i.gtf
	/home/emmanuel/software/gffcompare-0.11.5.Linux_x86_64/gffcompare -o ${folder}/$i/$i\_RQ -r ${ref_annot} -R -Q -D $i.gtf
done
