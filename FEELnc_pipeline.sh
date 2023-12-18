#!/bin/bash 

#conda activate STARalign

ref_gtf=$1
ref_fa=$2
gtf_list=$3
learning=$4

for i in `cat ${gtf_list}`
do
	mkdir FEELnc\_$i
	FEELnc_filter.pl -i $i.gtf -a ${ref_gtf} -b transcript_biotype=protein_coding > FEELnc\_$i/$i\_filtered.gtf
	FEELnc_codpot.pl -i FEELnc\_$i/$i\_filtered.gtf -a ${ref_gtf} -g ${ref_fa} -b transcript_biotype=protein_coding -n 20000,100 --outdir=FEELnc\_$i ${learning}
	FEELnc_classifier.pl -i FEELnc\_$i/$i.lncRNA.gtf -a ${ref_gtf} > FEELnc\_$i/$i\_classes.txt
done 


