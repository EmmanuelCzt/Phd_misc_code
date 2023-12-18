#!/bin/bash
# Alignement for downstream assembly using Scallop or Stringtie
# conda activate STARalign

sample_list=$1
idxgenome=$2
assembly=$3

mkdir map

# 1st Pass
for i in `cat ${sample_list}` 
do 
	STAR --runMode alignReads \
	--outSAMstrandField intronMotif \
	--outFilterType BySJout \
	--outFilterIntronMotifs RemoveNoncanonical \
	--alignEndsType EndToEnd \
	--genomeLoad LoadAndKeep \
	--readFilesCommand zcat \
	--outSAMtype BAM Unsorted \
	--genomeDir ${idxgenome} \
	--readFilesIn $i\_1.fastq.gz $i\_2.fastq.gz \
	--runThreadN 12 \
	--outFileNamePrefix map/$i\_${assembly}
done

#2nd Pass
#for i in `cat ${sample_list}` 
#do 
#	STAR --runMode alignReads \
#	--outSAMstrandField intronMotif \
#	--outFilterType BySJout \
#	--outFilterIntronMotifs RemoveNoncanonical \
#	--alignEndsType EndToEnd \
#	--genomeLoad NoSharedMemory \
#	--readFilesCommand zcat \
#	--outSAMtype BAM Unsorted \
#	--genomeDir ${idxgenome} \
#	--readFilesIn $i\_1.fastq.gz $i\_2.fastq.gz \
#	--runThreadN 12 \
#	--outFileNamePrefix map/$i\_P2\_${assembly} \
#	--sjdbFileChrStartEnd map/$i\_P1\_${assembly}SJ.out.tab
#done


STAR --genomeDir ${idxgenome} --genomeLoad Remove

# --outSAMstrandField intronMotif for compatibility with Stringtie and Scallop for unstranded datasets (XS information)
# --outFilterType BySJout  : filters out spurious junctions
# --outFilterIntronMotifs RemoveNoncanonical remove unnanotated non canonical junctions