#!/bin/bash

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
   hisat2 -p 8 --dta -x ${idxgenome} -1 $i/$i\_1.fastq.gz -2 $i/$i\_2.fastq.gz -S map/$i.sam |& tee -a paired_alignement.log
 done

cd map/

echo "***paired alignement done***" |& tee -a paired_alignement.log

# Sorting uniquely mapping reads for transcript reconstruction and 
echo "***generation of uniquely mapping reads files***" |& tee -a paired_alignement.log
for i in `cat ../${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools view -@ 8 -q 10 -F 1280 -h -b $i.sam > $i.F1280.bam |& tee -a paired_alignement.log # Remove multi mapping reads (flag : 256) and PCR duplicate (1024)
  samtools view -@ 8 -q 10 -F 1024 -h -b $i.sam > $i.F1024.bam |& tee -a paired_alignement.log # Impact of multi mapping reads on transcript reconstruction
done

echo "***bam generation done & deleting sam files" |& tee -a paired_alignement.log
rm *.sam

#sorting
echo "***sam sorting & bam generation***" 2>&1 | tee > paired_alignement.log
for i in `cat ../${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools sort -@ 8 -h -o $i.F1280.sorted.bam -O bam $i.F1280.bam |& tee -a paired_alignement.log
  samtools sort -@ 8 -h -o $i.F1024.sorted.bam -O bam $i.F1024.bam |& tee -a paired_alignement.log
done

# Merging bam files 
echo "***Merging***" |& tee -a paired_alignement.log
ls *.F1280.sorted.bam > F1280.list
ls *.F1024.sorted.bam > F1024.list

samtools merge -h -b F1280.list ${mergedfile}\_F1280.bam
samtools merge -h -b F1024.list ${mergedfile}\_F1024.bam

# Indexing and stats 
echo "***generation of stats files***" |& tee -a paired_alignement.log

samtools index -@ 8 -b ${mergedfile}\_F1280.bam |& tee -a paired_alignement.log
samtools index -@ 8 -b ${mergedfile}\_F1024.bam |& tee -a paired_alignement.log


samtools idxstats -@ 8 ${mergedfile}\_F1280.bam > F1280.stats |& tee -a paired_alignement.log
samtools flagstat -@ 8 ${mergedfile}\_F1280.bam > F1280.flagstat |& tee -a paired_alignement.log
samtools idxstats -@ 8 ${mergedfile}\_F1024.bam > F1024.stats |& tee -a paired_alignement.log
samtools flagstat -@ 8 ${mergedfile}\_F1024.bam > F1024.flagstat |& tee -a paired_alignement.log

cd ../

# Annotation
echo "***Annotation uising StringTie & Scallop***" |& tee -a paired_alignement.log
mkdir annotation

scallop -i map/${mergedfile}\_F1280.bam --min_transcript_length_base 200 --min_flank_length 10 --min_splice_bundary_hits 3 --min_bundle_gap 10 --min_transcript_coverage 0.1 --min_single_exon_coverage 2 -o annotation/${mergedfile}\_F1280.scallop.gtf |& tee -a paired_alignement.log
scallop -i map/${mergedfile}\_F1024.bam --min_transcript_length_base 200 --min_flank_length 10 --min_splice_bundary_hits 3 --min_bundle_gap 10 --min_transcript_coverage 0.1 --min_single_exon_coverage 2 -o annotation/${mergedfile}\_F1024.scallop.gtf |& tee -a paired_alignement.log

stringtie map/${mergedfile}\_F1280.bam -p 8 -f 0.50 -m 200 -a 10 -j 3 -c 0.1 -g 10 -o annotation/${mergedfile}\_F1280.stringtie.gtf |& tee -a paired_alignement.log
stringtie map/${mergedfile}\_F1024.bam -p 8 -f 0.50 -m 200 -a 10 -j 3 -c 0.1 -g 10 -o annotation/${mergedfile}\_F1024.stringtie.gtf |& tee -a paired_alignement.log


# generating bigWig files centered on the X chr
echo "***generation of X chromosome centered bigWigs***" |& tee -a paired_alignement.log
mkdir bigwig

bamCoverage --filterRNAstrand forward -r chrX --normalizeUsing BPM -b map/${mergedfile}\_F1280.bam -o bigwig/${mergedfile}\_F1280.forward.chrX.BPM.bw |& tee -a paired_alignement.log
bamCoverage --filterRNAstrand reverse -r chrX --normalizeUsing BPM -b map/${mergedfile}\_F1280.bam -o bigwig/${mergedfile}\_F1280.reverse.chrX.BPM.bw |& tee -a paired_alignement.log

bamCoverage --filterRNAstrand forward -r chrX --normalizeUsing BPM -b map/${mergedfile}\_F1024.bam -o bigwig/${mergedfile}\_F1024.forward.chrX.BPM.bw |& tee -a paired_alignement.log
bamCoverage --filterRNAstrand reverse -r chrX --normalizeUsing BPM -b map/${mergedfile}\_F1024.bam -o bigwig/${mergedfile}\_F1024.reverse.chrX.BPM.bw |& tee -a paired_alignement.log