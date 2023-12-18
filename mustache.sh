#!/bin/bash

INPUT=$1
CHROM=$2
RESOLUTION=$3
SPARSITY=$4
OUTDIR=$5

mkdir ${OUTDIR}
mkdir ${OUTDIR}/bedpe

#Call loops
mustache -f ${INPUT} -ch ${CHROM} -r ${RESOLUTION} -pt 0.05 -st ${SPARSITY} -o ${OUTDIR}/`basename ${INPUT} .mcool`_${RESOLUTION}.tsv

#make bedpe files
tail -n +2 ${OUTDIR}/`basename ${INPUT} .mcool`_${RESOLUTION}.tsv | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,"H9",$7,"+","+"}' > ${OUTDIR}/bedpe/`basename ${INPUT} .mcool`_${RESOLUTION}.bedpe

