#!/bin/bash

wigfile=$1
map_chain=$2

mkdir uniq

chmod +x /home/emmanuel/software/bigWigToBedGraph
chmod +x /home/emmanuel/software/wigToBigWig
chmod +x /home/emmanuel/software/bedGraphToBigWig
chmod +x /home/emmanuel/software/liftOver
chmod +x /home/emmanuel/software/bedRemoveOverlap

for i in `cat ${wigfile}`
do
  #/home/emmanuel/software/wigToBigWig $i.wig ${chrom_sizes_old} $i\_temp.bigWig
  /home/emmanuel/software/bigWigToBedGraph $i.bigWig $i.bedGraph
  grep "^chrX" $i.bedGraph | sort -k1,1 -k2,2n  > $i.chrX.sorted.bedGraph
  /home/emmanuel/software/bedRemoveOverlap $i.chrX.sorted.bedGraph $i.chrX.sorted.no_overlap.bedGraph
  /home/emmanuel/software/liftOver -bedPlus=4 -positions $i.chrX.sorted.no_overlap.bedGraph ${map_chain} $i.chrX.sorted.no_overlap.lifted.bedGraph $i.unmapped
  sort -k1,1 -k2,2n $i.chrX.sorted.no_overlap.lifted.bedGraph > $i.chrX.sorted2.no_overlap.lifted.bedGraph
  /home/emmanuel/software/bedRemoveOverlap $i.chrX.sorted2.no_overlap.lifted.bedGraph $i.chrX.sorted2.no_overlap2.lifted.bedGraph
  cut -f 1 | sort $i.chrX.sorted2.no_overlap2.lifted.bedGraph > $i.col1
  comm -23 $i.col1 rheMac8.sorted.col1  > uniq/$i\_col1.uniq
done
