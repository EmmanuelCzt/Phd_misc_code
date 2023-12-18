#Get only assembled chromosomes
grep -Ew "^chr([1-9]|[1-9][0-9]|[X-Y])"

#Produce exon bed file from GTF and remove the duplicates
awk -v OFS="\t" '$3=="exon" {print $1,$4,$5,$10,".",$7}'  Macaca_mulatta.Mmul_10.108.UCSC.XICRNA.20221129.gtf | tr -d '";' | sort -k1,1 -k2,2n | uniq -u > Macaca_mulatta.Mmul_10.108.UCSC.XICRNA.20221129.exons.unique.bed

sort -k1,1V -k2,2n -k3,3n