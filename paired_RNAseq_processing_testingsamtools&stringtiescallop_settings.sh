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
mergedfile=$3

mkdir map

echo "files list : `cat ${sample_list}`" |& tee -a paired_alignement.log
echo "genome used : ${idxgenome}" |& tee -a paired_alignement.log

echo "***paired alignement using hisat2***" |& tee -a paired_alignement.log
# Paired alignement with 16 threads and dta for StringTie
for i in `cat ${sample_list}`
do
   echo $i
   hisat2 -p 8 --dta -x ${idxgenome} -1 $i\_R1.fastq.gz -2 $i\_R2.fastq.gz -S map/$i.sam |& tee -a paired_alignement.log
 done

echo "***paired alignement done***" |& tee -a paired_alignement.log

# Sorting uniquely mapping reads for transcript reconstruction and 
echo "***generation of uniquely mapping reads files & removal of multimapping/PCR duplicates reads***" |& tee -a paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools view -@ 8 -q 10 -F 1280 -h -b map/$i.sam > map/$i.F1280.bam |& tee -a paired_alignement.log # Remove multi mapping reads (flag : 256) and PCR duplciate (1024)
done

echo "***bam generation done & deleting sam files" |& tee -a paired_alignement.log
rm -r map/*.sam

#sort and formatting to bam & 8 threads
echo "***bam sorting***" 2>&1 | tee > paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools sort -@ 8 -h -o map/$i.F1280.sorted.bam map/$i.F1280.bam |& tee -a paired_alignement.log
done

# generating stat files : Pre clean-up
echo "***generation of stats files***" |& tee -a paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools idxstats -@ 8 map/$i.sorted.bam > $i.stats |& tee -a paired_alignement.log
  samtools flagstat -@ 8 map/$i.sorted.bam > $i.flagstat |& tee -a paired_alignement.log
done

# generating index files, useful for IGV visualization

echo "***generation of index files***" |& tee -a paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools index -@ 8 -b map/$i.sorted.bam |& tee -a paired_alignement.log
done

# Sorting uniquely mapping reads for transcript reconstruction and 
echo "***generation of uniquely mapping reads files***" |& tee -a paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools view -@ 8 -q 10 -F 1280 -h -b map/$i.sorted.bam > map/$i.sorted.F1280.bam |& tee -a paired_alignement.log # Remove multi mapping reads (flag : 256) and PCR duplciate (1024)
  samtools view -@ 8 -q 10 -F 1024 -h -b map/$i.sorted.bam > map/$i.sorted.F1024.bam |& tee -a paired_alignement.log # Impact of multi mapping reads on transcript reconstruction
done

# generating stat files : Pre clean-up
echo "***generation of stats files***" |& tee -a paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools idxstats map/$i.sorted.bam > $i.clean.stats |& tee -a paired_alignement.log
  samtools flagstat map/$i.sorted.bam > $i.clean.flagstat |& tee -a paired_alignement.log
done

cd map/

# Merging .bam files
echo "***Merging***" |& tee -a paired_alignement.log
ls *.sorted.F1280.bam > F1280.list
ls *.$i.sorted.F1024.bam > F1024.list

samtools merge -h -b F1280.list ${mergedfile}\_F1280.bam |& tee -a paired_alignement.log
samtools merge -h -b F1024.list ${mergedfile}\_F1024.bam |& tee -a paired_alignement.log

# Sorting, indexing and stats
sammtools sort -@ 8 -h -b 


# Annotation
echo "***Annotation uising StringTie & Scallop***" |& tee -a paired_alignement.log
mkdir annotation

for i in `cat ${sample_list}`
do
  scallop -i map/$i.sorted.bam -o annotation/$i.scallop.gtf |& tee -a paired_alignement.log
  stringtie map/$i.sorted.bam -p 8 -o annotation/$i.stringtie.gtf |& tee -a paired_alignement.log
done
# generating bigWig files centered on the X chr
echo "***generation of X chromosome centered bigWigs***" |& tee -a paired_alignement.log

mkdir bigwig

for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  bamCoverage -p 8 --filterRNAstrand forward -r chrX --normalizeUsing BPM -b map/$i.sorted.bam -o bigwig/$i.forward.chrX.BPM.bw |& tee -a paired_alignement.log
  bamCoverage -p 8 --filterRNAstrand reverse -r chrX --normalizeUsing BPM -b map/$i.sorted.bam -o bigwig/$i.reverse.chrX.BPM.bw |& tee -a paired_alignement.log
done
