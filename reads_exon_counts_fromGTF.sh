#!/bin/bash

sample_list=$1

echo "Formatting gtf to bed with count values"
for i in `cat ${sample_list}`
do
	
	awk '$3=="transcript" {print}' $i.gtf | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$10,$12,$14}' | sed -E 's/"|;//g' | sed 's/ /\t/g' >count=0; second=0;for i in `cat chr.list`; do count=$((count+1)) ; second=$((${count}+1)) ; echo $second ; $i.FPKM
done

echo "counting exons per transcript"
for i in `cat ${sample_list}`
do
	for j in `awk '{print $10}' $i.FPKM | grep -v "transcript_id"` # Using count file is better since transcript_ids are already formatted 
	do 
		
		awk '{print $12}' $i.gtf | grep -Fo "$j" | tail -n +2 | awk 'END {print NR}' >> $i.exon # grep -F to take input as literal pattern since there are "." in the ranscript ids

	done
	paste $i.FPKM $i.exon > $i.extra #tsv file with extra information on reads supporting each transcript and its number of exons 

done
 	
#tsv file
#chr	method	feature	start	end score strand smthg	gene_id transcript_id 	rpkm exon_counts