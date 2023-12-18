#!/bin/bash

wigfile=$1
chrom_sizes_old=$2
chrom_sizes_new=$3
map_chain=$4

chmod +x /home/emmanuel/software/bigWigToBedGraph
chmod +x /home/emmanuel/software/wigToBigWig
chmod +x /home/emmanuel/software/bedGraphToBigWig
chmod +x /home/emmanuel/software/liftOver
chmod +x /home/emmanuel/software/bedRemoveOverlap

for i in `cat ${wigfile}`
do
  /home/emmanuel/software/wigToBigWig $i.wig ${chrom_sizes_old} $i\_temp.bw
  /home/emmanuel/software/bigWigToBedGraph $i\_temp.bw $i.bedGraph
  grep "^chrX" $i.bedGraph | sort -k1,1 -k2,2n  > $i.chrX.sorted.bedGraph
  /home/emmanuel/software/bedRemoveOverlap $i.chrX.sorted.bedGraph $i.chrX.sorted.no_overlap.bedGraph
  /home/emmanuel/software/liftOver -bedPlus=4 -positions $i.chrX.sorted.no_overlap.bedGraph ${map_chain} $i.chrX.sorted.no_overlap.hg38.bedGraph $i.unmapped
  sort -k1,1 -k2,2n $i.chrX.sorted.no_overlap.hg38.bedGraph > $i.chrX.sorted2.no_overlap.hg38.bedGraph
  /home/emmanuel/software/bedRemoveOverlap $i.chrX.sorted2.no_overlap.hg38.bedGraph $i.chrX.sorted2.no_overlap2.hg38.bedGraph
  /home/emmanuel/software/bedGraphToBigWig $i.chrX.sorted2.no_overlap2.hg38.bedGraph ${chrom_sizes_new} $i.chrX.bw
done

rm *.bedGraph
rm *.unmapped
rm *temp.bw
