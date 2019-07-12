export PATH=$PATH:/home/emmanuel/software/sratoolkit.2.9.6-1-ubuntu64/bin
export PATH=$PATH:/home/user/software/hisat2-2.1.0-Linux_x86_64/hisat2-2.1.0 # reads alignement
export PATH=$PATH:/home/user/software/samtools # bam generation and sorting
export PATH=$PATH:/home/user/software/scallop-0.10.3_linux_x86_64 # transcript annotation
export PATH=$PATH:/home/user/software/stringtie-1.3.5.Linux_x86_64 # transcript annotation

sample_list=$1

for i in `cat ${sample_list}`
do
  mkdir $i
  fasterq-dump $i -O $i/
  gzip $i/$i\_1.fastq
  gzip $i/$i\_2.fastq
  rm -r $i/*.fastq
done

idxgenome=$2 # genome assembly index generated using : hisat2-build genome.fa output_filename

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

echo "***paired alignement done***" |& tee -a paired_alignement.log

# sort and formatting to bam & 8 threads
echo "***sam sorting & bam generation***" 2>&1 | tee > paired_alignement.log
for i in `cat ${sample_list}`
do
  echo $i |& tee -a paired_alignement.log
  samtools sort -@ 8 -o map/$i.sorted.bam map/$i.sam |& tee -a paired_alignement.log
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
