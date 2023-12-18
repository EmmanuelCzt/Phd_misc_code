#!/bin/bash

sample_list=$1
read_size=$2


mkdir trimm

for i in `cat ${sample_list}`
do
  echo $i | tee > cutadapt.log
  cutadapt -q 20 -m ${read_size} -b AAGCAGTGGTATCAACGCAGAGTACT -B AAGCAGTGGTATCAACGCAGAGTACT -a T{50} -A T{50} -j 8 -o trimm/$i\_1.trimmed.fastq.gz -p trimm/$i\_2.trimmed.fastq.gz $i\_1.fastq.gz $i\_2.fastq.gz 2>&1 | tee > cutadapt.log
done

# -q : quality score base 33
# -m : minimum length of the trimmed reads
# -j : multi threading (pigz package needed)
