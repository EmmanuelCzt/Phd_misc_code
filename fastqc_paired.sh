#!/bin/bash

#export PATH=$PATH:/home/emmanuel/software/fastqc_v0.11.8/FastQC

sample_list=$1

mkdir QC
mkdir QC/multiqc

for i in `cat ${sample_list}`
do
  #echo $i 2>&1 | tee > fastqc/paired_alignement.log
  fastqc -t 8 -o QC/ fastq/$i\_1.fastq.gz fastq/$i\_2.fastq.gz
done

multiqc -o QC/multiqc QC/