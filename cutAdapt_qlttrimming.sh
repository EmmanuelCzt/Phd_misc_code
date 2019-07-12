#!/bin/bash

sample_list=$1

export PATH=$PATH:/home/emmanuel/software/fastqc_v0.11.8/FastQC

mkdir trimmed

for i in `cat ${sample_list}`
do
  echo $i | tee > cutadapt.log
  cutadapt -q 30 -m 50 -j 8 -o trimmed/$i\_1.trimmed.fastq.gz -p trimmed/$i\_2.trimmed.fastq.gz $i\_1.fastq.gz $i\_2.fastq.gz 2>&1 | tee > cutadapt.log
  fastqc -o fastqc/ rimmed/$i\_1.trimmed.fastq.gz trimmed/$i\_2.trimmed.fastq.gz 2>&1 | tee > cutadapt.log
done

# -q : quality score base 33
# -m : minimum length of the trimmed reads
# -j : multi threading (pigz package needed)
