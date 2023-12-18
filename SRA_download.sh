#!/bin/bash

#export PATH=$PATH:/home/emmanuel/software/sratoolkit.current-ubuntu64/sratoolkit.2.9.6-1-ubuntu64/bin/

sra_list=$1

mkdir fastq

for i in `cat ${sra_list}`
do
  
  prefetch -p -O . $i
  fasterq-dump $i.sra -O fastq/
done

gzip fastq/*.fastq

rm *.sra
