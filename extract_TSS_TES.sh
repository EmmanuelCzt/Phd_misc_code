#!/bin/bash

GTF=$1
OUTPREFIX=$2

mkdir ${OUTPREFIX}

##Extract TSS for sense and antisense transcripts

#Transcript id column is 12 for gencode, 14 for ensemble gtf files

awk '$3=="transcript" {print $0}' ${GTF}| awk -v OFS="\t" '$7=="+" {print $1,$4,$4+1,$14,0,$7}' | \
tr -d '";' > ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TSS_fwd.bed

awk '$3=="transcript" {print $0}' ${GTF} | awk -v OFS="\t" '$7=="-" {print $1, $5,$5-1,$14,0,$7}' | \
tr -d '";' > ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TSS_rev.bed

#Merge and sort forward and reverse TSS
cat ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TSS_fwd.bed ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TSS_rev.bed | \
sort -k1,1 -k2,2n > ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TSS.bed


##Extract TES for sense and antisense transcripts
awk '$3=="transcript" {print $0}' ${GTF} | awk -v OFS="\t" '$7=="+" {print $1,$5-1,$5,$14,0,$7}' | \
tr -d '";' > ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TES_fwd.bed

awk '$3=="transcript" {print $0}' ${GTF} | awk -v OFS="\t" '$7=="-" {print $1,$4,$4-1,$14,0,$7}' | \
tr -d '";' > ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TES_rev.bed

#Merge and sort forward and reverse TES
cat ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TES_fwd.bed ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TES_rev.bed | \
sort -k1,1 -k2,2n > ${OUTPREFIX}/${OUTPREFIX}\_transcripts_TES.bed