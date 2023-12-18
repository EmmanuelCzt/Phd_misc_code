#!/bin/bash

REF=$1 # coordinates reference File to be mapped to
DATA=$2 # cordinates file with mesure for a given genomic interval
WINSIZE=$3 #window sizes 
OUTNAME=$4

#Mapping sequencing signal around features

#Make windows
bedtools makewindows -b ${REF} -w ${WINSIZE} -i winnum \
| sort -k1,1 -k2,2n | tr "_" "\t" > `basename ${REF} .bed`.${WINSIZE}bpwindows.bed

#map data to ref
bedtools map -a `basename ${REF} .bed`.${WINSIZE}bpwindows.bed \
-b ${DATA} -c 4 -o collapse -null 0 > ${OUTNAME}.map.bed

#Summerize mapped data by windows
sort -t$'\t' -k5,5n ${OUTNAME}.map.bed | \
bedtools groupby -i - -g 4 -c 5 -o collapse -null 0 > ${OUTNAME}.${WINSIZE}bp.counts.txt