#!/bin/bash

list=$1

for i in `cat ${list}`
do
	awk '$3=="transcript" {print $14}' $i\_shuffle/$i\.gtf.lncRNA.gtf | tr -d '";' | uniq | sort > $i\_shuffle/$i\_lncRNA_sorted.txt
	awk '$3=="transcript" {print $14}' $i\_shuffle/$i\.gtf.mRNA.gtf | tr -d '";' | uniq | sort > $i\_shuffle/$i\_mRNA_sorted.txt
	comm -12 $i\_shuffle/$i\_lncRNA_sorted.txt $i\_shuffle/$i\_mRNA_sorted.txt > $i\_shuffle/$i\_shared.txt
	comm -23 $i\_shuffle/$i\_lncRNA_sorted.txt $i\_shuffle/$i\_mRNA_sorted.txt > $i\_shuffle/$i\_lncRNA.txt
	comm -13 $i\_shuffle/$i\_lncRNA_sorted.txt $i\_shuffle/$i\_mRNA_sorted.txt > $i\_shuffle/$i\_mRNA.txt
done 