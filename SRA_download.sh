#!/bin/bash

export PATH=$PATH:/home/emmanuel/software/sratoolkit.2.9.6-ubuntu64/bin

sra_list=$1

for i in `cat ${sra_list}`
do
  mkdir $i
  fasterq-dump $i -O $i/
  gzip $i/$i\_1.fastq
  gzip $i/$i\_2.fastq
  rm -r $i/*.fastq 
done
