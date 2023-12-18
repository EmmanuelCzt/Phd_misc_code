#!/bin/bash

export SpeciesTree=$1

export path=/mnt/mydatalocal/IPLOSS
export pathGenomes=${path}/data/genome_sequences
export pathResults=${path}/results/whole_genome_alignments/

mkdir -p ${pathResults}/tmp

#########################################################################
# Create params file
cp ${SpeciesTree} ${pathResults}/seqFile # Species tree

echo "" >>  ${pathResults}/seqFile

# Species list and path to genome file
for sp in `ls ${pathGenomes}`
do
    echo "${sp} ${pathGenomes}/${sp}/genome_sm.fa" >> ${pathResults}/seqFile
done

#########################################################################

docker run -v ${path}:/mnt/mydatalocal/IPLOSS --rm -t quay.io/comparative-genomics-toolkit/cactus:v2.4.0 cactus --maxCores 32 --maxMemory 100G --maxDisk 500G --defaultDisk 20G --defaultMemory 4G --binariesMode local --workDir ${pathResults}/ ${pathResults}/jobStore ${pathResults}/seqFile ${pathResults}/alignment.hal

#########################################################################
