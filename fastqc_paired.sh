#!/bin/bash

export PATH=$PATH:/home/emmanuel/software/fastqc_v0.11.8/FastQC

sample_list=$1

mkdir fastqc

for i in `cat ${sample_list}`
do
  echo $i 2>&1 | tee > fastqc/paired_alignement.log
  fastqc -o fastqc/ $i\_1.fastq.bz2 $i\_2.fastq.bz2 2>&1 | tee > fastqc/paired_alignement.log
done
