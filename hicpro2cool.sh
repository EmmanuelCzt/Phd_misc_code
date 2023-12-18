#!/bin/bash

##Convert HiCpro files to cool format
##Permorf ice normalization
##plot un balancer and balanced conact frequencies


sample_list=$1
bin_sizes=$2
chrom_sizes=$3
region=$4

mkdir cool
mkdir cool_pics


for j in `cat ${bin_sizes}`
do
	cooler load -f coo --one-based ${chrom_sizes}:$j raw/${sample_list}\_$j.matrix cool/${sample_list}\_$j.cool
	cooler balance -p 10 cool/${sample_list}\_$j.cool
	cooler show --out cool_pics/${sample_list}\_$j\_cool.svg --dpi 200 cool/${sample_list}\_$j.cool ${region}
	cooler show -b --out cool_pics/${sample_list}\_$j\_cool_balanced.svg --dpi 200 cool/${sample_list}\_$j.cool ${region}
done