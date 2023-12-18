#!/bin/bash
# Alignement for downstream assembly using Scallop or Stringtie
# conda activate STARalign
sj=$1
	STAR --runMode alignReads \
	--outSAMstrandField intronMotif \
	--outFilterType BySJout \
	--outFilterIntronMotifs RemoveNoncanonical \
	--alignEndsType EndToEnd \
	--genomeLoad NoSharedMemory \
	--readFilesCommand zcat \
	--outSAMtype BAM Unsorted \
	--genomeDir /home/emmanuel/Documents/Genomes/hg38/STARindex_hg38_Ens99 \
	--readFilesIn SRR3192400_1.fastq.gz,SRR3192402_1.fastq.gz,SRR3192414_1.fastq.gz,SRR3192416_1.fastq.gz,SRR3192544_1.fastq.gz,SRR3192544_1.fastq.gz,SRR4422432_1.fastq.gz SRR3192400_2.fastq.gz,SRR3192402_2.fastq.gz,SRR3192414_2.fastq.gz,SRR3192416_2.fastq.gz,SRR3192544_2.fastq.gz,SRR3192544_2.fastq.gz,SRR4422432_2.fastq.gz \
	--runThreadN 12 \
	--outFileNamePrefix map/P2 \
	--sjdbFileChrStartEnd `cat ${sj}`