#!/bin/bash
# Align stranded paired RNA-seq reads with downstream applications as inference of a transcript model from the data and coverage plots.
# Usage : ./paired_RNAseq_processing sample_name hisat_genome_idx
#Add software to PATH
#export PATH=$PATH:/home/user/software/hisat2-2.1.0-Linux_x86_64/hisat2-2.1.0 # reads alignement
#export PATH=$PATH:/home/user/software/samtools # bam generation and sorting
#export PATH=$PATH:/home/user/software/scallop-0.10.3_linux_x86_64 # transcript annotation
#export PATH=$PATH:/home/user/software/stringtie-1.3.5.Linux_x86_64 # transcript annotation

# Variables : fastq (if multiple files make a list) idxgenome
sample_list=$1 
mkdir map

echo "***paired alignement done***" 2>&1 >> samtools_ops.log

# Sorting uniquely mapping reads (q<10) & filtering multi mapping reads (flag : 256) & PCR duplicates PCR duplciate (1024)
echo "***generation of uniquely mapping reads files & removal of multimapping/PCR duplicates reads***" 2>&1 >> samtools_ops.log
samtools --version 2>&1 >> samtools_ops.log

for i in `cat ${sample_list}`
do
  echo $i 2>&1 >> samtools_ops.log
  samtools view -@ 8 -q 10 -F 1280 -h -b map/$i.sam > map/$i.F1280.bam
done

echo "***bam generation done & deleting sam files" 2>&1 >> samtools_ops.log
rm -r map/*.sam

#sort and formatting to bam & 8 threads
echo "***bam sorting***"  2>&1 > samtools_ops.log
for i in `cat ${sample_list}`
do
  echo $i 2>&1 >> samtools_ops.log
  samtools sort -@ 8 -o map/$i.F1280.sorted.bam map/$i.F1280.bam 2>&1 >> samtools_ops.log
done

rm -r map/*.F1280.bam

# generating index files, useful for IGV visualization
echo "***generation of index files***" 2>&1 >> samtools_ops.log
for i in `cat ${sample_list}`
do
  echo $i 2>&1 >> samtools_ops.log
  samtools index -@ 8 -b map/$i.F1280.sorted.bam 2>&1 >> samtools_ops.log
done

# generating stat files

mkdir map/stats

echo "***generation of stats files***" 2>&1 >> samtools_ops.log
for i in `cat ${sample_list}`
do
  echo $i 2>&1 >> samtools_ops.log
  samtools idxstats map/$i.F1280.sorted.bam > map/stats/$i.stats 2> samtools_ops.log
  samtools flagstat map/$i.F1280.sorted.bam > map/stats/$i.flagstat 2> samtools_ops.log
done


