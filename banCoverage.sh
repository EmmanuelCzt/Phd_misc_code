#!/bin/bash

sample_list=$1

 mkdir bigwig

 for i in `cat ${sample_list}`
 do
  echo $i | tee >> paired_alignement.log
  bamCoverage -p 8 --filterRNAstrand forward -r  chrX:65178685:66401555 --normalizeUsing BPM -b map/$i.sorted.bam -o bigwig/$i.forward.chrX.BPM.bw 2>&1 | tee >> paired_alignement.log
  bamCoverage -p 8 --filterRNAstrand reverse -r  chrX:65178685:66401555 --normalizeUsing BPM -b map/$i.sorted.bam -o bigwig/$i.reverse.chrX.BPM.bw 2>&1 | tee >> paired_alignement.log
 done
