#!/bin/bash

#List of Gene ID to grep
GENEID=$1

#empirical Tx model from Scallop, TACO or stringtie
EMPTX=$2
#Reference transcript model to add the Tx to 
REFTX=$3

#outfile 
OUTFILE=$4

mkdir temp

for i in `cat ${GENEID}`
do
	grep $i ${EMPTX} >> temp/subset_Tx.gtf
done

#merge transcripts
stringtie --merge -G ${REFTX} -o temp/stringtie.merge.gtf temp/subset_Tx.gtf

#Regenerate full ref transcript model
awk -v OFS="\t" '$3!="transcript" && $3!="exon" {print $0}' ${REFTX} > temp/features.gtf

cat temp/stringtie.merge.gtf temp/features.gtf > ${OUTFILE}

rm -r temp/