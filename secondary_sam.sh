#!/bin/bash

input=$1

mkdir scd_almt

samtools view -f 256 ${input} | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$12,$13}' | sed -e 's/[NM:i:|MD:Z:]//g' > scd_almt/`basename ${input} .bam`_f256.samout
samtools view -F 256 ${input} | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$12,$13}' | sed -e 's/[NM:i:|MD:Z:]//g' > scd_almt/`basename ${input} .bam`_F256.samout