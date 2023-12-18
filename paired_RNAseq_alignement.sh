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
idxgenome=$2 # genome assembly index generated using : hisat2-build genome.fa output_filename

mkdir map

echo "file list : `cat ${sample_list}`" 2>&1 >> paired_alignement.log
echo "genome used : ${idxgenome}" 2>&1 >> paired_alignement.log

echo "***paired alignement using hisat2***" 2>&1 >> paired_alignement.log
hisat2 --version 2>&1 >> paired_alignement.log
# Paired alignement with 16 threads and dta for StringTie
for i in `cat ${sample_list}`
do
   echo $i
   hisat2 -p 8 --dta -x ${idxgenome} -1 fastq/$i\_R1.fastq.gz -2 fastq/$i\_R2.fastq.gz -S map/$i.sam 2>&1 >> paired_alignement.log
 done

echo "***paired alignement done***" 2>&1 >> paired_alignement.log

# Sorting uniquely mapping reads (q<10) & filtering multi mapping reads (flag : 256) & PCR duplicates PCR duplciate (1024)
echo "***generation of uniquely mapping reads files & removal of multimapping/PCR duplicates reads***" 2>&1 >> paired_alignement.log
samtools --version 2>&1 >> paired_alignement.log

for i in `cat ${sample_list}`
do
  echo $i 2>&1 >> paired_alignement.log
  samtools view -@ 8 -q 10 -F 1280 -h -b map/$i.sam > map/$i.F1280.bam 2>> paired_alignement.log
done

echo "***bam generation done & deleting sam files" 2>&1 >> paired_alignement.log
rm -r map/*.sam

#sort and formatting to bam & 8 threads
echo "***bam sorting***" 2>&1 >> paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i 2>&1 >> paired_alignement.log
  samtools sort -@ 8 -o map/$i.F1280.sorted.bam map/$i.F1280.bam 2>&1 >> paired_alignement.log
done

# generating index files, useful for IGV visualization
echo "***generation of index files***" 2>&1 >> paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i 2>&1 >> paired_alignement.log
  samtools index -@ 8 -b map/$i.F1280.sorted.bam 2>&1 >> paired_alignement.log
done

# generating stat files

mkdir map/stats

echo "***generation of stats files***" 2>&1 >> paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i 2>&1 >> paired_alignement.log
  samtools idxstats map/$i.F1280.sorted.bam > map/stats/$i.stats 2>> paired_alignement.log
  samtools flagstat map/$i.F1280.sorted.bam > map/stats/$i.flagstat 2>> paired_alignement.log
done

# Annotation
echo "***Annotation uising StringTie & Scallop***" 2>> paired_alignement.log
mkdir annotation

for i in `cat ${sample_list}`
do
  scallop -i map/$i.F1280.sorted.bam --min_transcript_length_base 200 --min_flank_length 10 --min_splice_bundary_hits 3 --min_bundle_gap 10 --min_transcript_coverage 0.1 --min_single_exon_coverage 20 -o annotation/$i.F1280.sorted.${organisme}.scallop.gtf |& tee -a paired_alignement.log
done
# generating bigWig files centered on the X chr
echo "***generation of X chromosome centered bigWigs***" 2>> paired_alignement.log

mkdir bigwig

for i in `cat ${sample_list}`
do
  echo $i 2>> paired_alignement.log
  bamCoverage -p 8 --filterRNAstrand forward -r chrX --normalizeUsing BPM -b map/$i.F1280.sorted.bam -o bigwig/$i.forward.chrX.BPM.bw 2>> paired_alignement.log
  bamCoverage -p 8 --filterRNAstrand reverse -r chrX --normalizeUsing BPM -b map/$i.F1280.sorted.bam -o bigwig/$i.reverse.chrX.BPM.bw 2>> paired_alignement.log
done

