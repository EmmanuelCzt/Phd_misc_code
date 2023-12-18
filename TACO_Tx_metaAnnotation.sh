#!/bin/bash


GTF=$1 #List of GTF
OUTDIR=$2

TACO=("/home/emmanuel/software/taco-v0.7.3.Linux_x86_64/taco_run")

${TACO} -v -o ${OUTDIR} --gtf-expr-attr cov --filter-min-expr 1.0 --isoform-frac 0.05 ${GTF}
