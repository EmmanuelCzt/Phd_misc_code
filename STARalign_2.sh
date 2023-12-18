#!/bin/bash
# If paired-end RNA-SEQ fastq1 list : coma separated file names. Fastq2 list : space separated file names
# conda activate STARalign

sample_list=$1
idxgenome=$2
assembly=$3
comp=$4
#SJ_tab=$5 # two pass mapping for SNP calling


mkdir map


STAR --runMode alignReads \
--outSAMstrandField intronMotif \
--outFilterType BySJout \
--outFilterIntronMotifs RemoveNoncanonical \
--readFilesCommand ${comp} \
--outSAMtype BAM Unsorted \
--genomeDir ${idxgenome} \
--readFilesIn `cat ${sample_list}` \
--runThreadN 12 \
--outFileNamePrefix map/${assembly} \
	



#PS RNA-seq data with a lot of unmapped reads so I will keep them to see what's up
#--outReadsUnmapped Fastx \
#--sjdbFileChrStartEnd `cat ${SJ_tab}`

