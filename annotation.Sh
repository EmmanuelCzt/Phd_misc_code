#!/bin/bash

# Strintie and scallop to path
export PATH=$PATH:/home/emmanuel/software/scallop-0.10.3_linux_x86_64
export PATH=$PATH:/home/emmanuel/software/stringtie-1.3.5.Linux_x86_64

sample_list=$1

mkdir annotation

for i in `cat ${sample_list}`
do
  scallop -i map/$i.sorted.bam -o annotation/$i.scallop.gtf
  stringtie map/$i.sorted.bam -p 4 -l stringtie -o annotation/$i.stringtie.gtf
done
