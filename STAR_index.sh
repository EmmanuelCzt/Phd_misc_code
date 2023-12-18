#!/bin/bash

outDir=$1
genomeFa=$2
# gtf=$3 --sjdbGTFfile ${gtf}

mkdir ${outDir}

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ${outDir} --readFilesCommand zcat --genomeFastaFiles ${genomeFa}