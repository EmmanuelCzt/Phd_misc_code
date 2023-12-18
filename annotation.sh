#!/bin/bash

mkdir annotation

sample_list=$1
organisme=$2

for i in `cat ${sample_list}`
do 
	scallop -i map/$i.F1280.sorted.bam --min_transcript_length_base 200 --min_flank_length 10 --min_splice_bundary_hits 3 --min_bundle_gap 50 --min_transcript_coverage 0.1 --min_single_exon_coverage 20 -o annotation/$i.F1280.sorted.${organisme}.scallop.gtf
	stringtie map/$i.F1280.sorted.bam -p 8 -f 0.01 -m 200 -a 10 -j 3 -c 0.1 -g 50 -o annotation/$i.F1280.sorted.${organisme}.stringtie.gtf
done

