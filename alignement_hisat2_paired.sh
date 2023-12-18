#!/bin/bash
# Align stranded paired RNA-seq reads with downstream applications as inference of a transcript model from the data and coverage plots.
# Usage : ./paired_RNAseq_processing sample_name hisat_genome_idx
#Add software to PATH
export PATH=$PATH:/home/emmanuel/software/hisat2-2.1.0-Linux_x86_64/hisat2-2.1.0 # reads alignement
export PATH=$PATH:/home/emmanuel/software/samtools # bam generation and sorting

# Variables : fastq (if multiple files make a list) idxgenome
sample_list=$1
idxgenome=$2 # genome assembly index generated using : hisat2-build genome.fa output_filename

mkdir map

echo "files list : `cat ${sample_list}`" |& tee -a paired_alignement.log
echo "genome used : ${idxgenome}" |& tee -a paired_alignement.log

echo "***paired alignement using hisat2***" |& tee -a paired_alignement.log
# Paired alignement with 16 threads and dta for StringTie
for i in `cat ${sample_list}`
do
   echo $i
   hisat2 -p 16 --dta -x ${idxgenome} -1 $i/$i\_1.fastq.gz -2 $i/$i\_2.fastq.gz -S map/$i.sam |& tee -a paired_alignement.log
 done

echo "***paired alignement done***" |& tee -a paired_alignement.log

# sort and formatting to bam & 8 threads
echo "***sam sorting & bam generation***" 2>&1 | tee > paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools sort -@ 16 -o map/$i.sorted.bam map/$i.sam |& tee -a paired_alignement.log
done

echo "***bam generation done & deleting sam files" |& tee -a paired_alignement.log
rm -r map/*.sam

# generating index files, useful for IGV visualization

echo "***generation of index files***" |& tee -a paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools index -b map/$i.sorted.bam |& tee -a paired_alignement.log
done

# Sorting uniquely mapping reads for counting
echo "***generation of uniquely mapping reads files***" |& tee -a paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools view -q10 -b map/$i.sorted.bam > map/$i.sorted.uniq.bam |& tee -a paired_alignement.log
done

# generating stat files (useful for samples sexing : readsChrX/readsChrY)
echo "***generation of stats files***" |& tee -a paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools idxstats map/$i.sorted.bam > $i.stats |& tee -a paired_alignement.log
done
