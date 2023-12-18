#!/bin/bash
# If paired-end RNA-SEQ fastq1 list : coma separated file names. Fastq2 list : space separated file names
# conda activate STARalign

sample=$1
idxgenome=$2
assembly=$3
comp=$4 #zcat of bzcat
outdir=$5
#SJ_tab=$5 # two pass mapping for SNP calling

mkdir ${outdir}
for i in `cat ${sample}`
do
	STAR --runMode alignReads \
	--genomeLoad LoadAndKeep \
	--outSAMstrandField intronMotif \
	--outFilterType BySJout \
	--outFilterIntronMotifs RemoveNoncanonical \
	--outSAMattributes NH HI AS nM NM \
	--readFilesCommand ${comp} \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir ${idxgenome} \
	--readFilesIn $i.fastq.gz \
	--runThreadN 6 \
	--limitBAMsortRAM 15000000000 \
	--outFileNamePrefix ${outdir}/`basename $i`\_${assembly} 
done


STAR --genomeDir ${idxgenome} --genomeLoad Remove

#PS RNA-seq data with a lot of unmapped reads so I will keep them to see what's up
#	--outReadsUnmapped Fastx \
#	--sjdbFileChrStartEnd ${SJ_tab}




