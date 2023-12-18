#!/usr/bin/env Rscript

library(ggplot2)
library(GenomicFeatures)
library(stringr)
setwd("/home/emmanuel/Documents/capture_HiC/SNP_id/mpileup/human/H9primedMC/")
#argmts <- commandArgs(trailingOnly = T) # chrom sizes or gtf = arg 2, boundaries = arg 3&4, DpnII file = arg 5

txdb <- makeTxDbFromGFF("/home/emmanuel/Documents/annotations/hg38/gencode.v32.annotation.gtf", format = "gtf")

#genomic boundaries
lbound <- c(72413461)
rbound <- c(75413461)

#Load DpnII file
#dpn2 <- read.table(argmts[4], sep = "\t", dec = ".",  header = F, col.names = c("id","chr","start","end","size","%GC"))

# Loading datasets, must be in working directory
col_id <- c("chr","pos","id","ref","alt","qual","filter","infos","genotype","GT-ll")
fileslistSNP <- list.files(".", pattern = "*[7|8|X].vcf")
filesSNPs <- lapply(fileslistSNP, function(x)read.table(x, sep = "\t", dec = ".", stringsAsFactors = F, header = F, 
                                                        col.names = col_id))
names(filesSNPs) <- c("chr7","chr8", "chrX")

# Extracting reads infos from each SNPs pos and filtering
extrac_DP4 <- function(x) { #Extract reads number for each SNP DP4=Ref+,Ref-,Alt+,Alt-
  unlist(unlist(str_extract_all(unlist(strsplit(x,split=";")),"DP4=[0-9]+,[0-9]+,[0-9]+,[0-9]+")))
}

extract_reads <- function(x) { #Extract read numbers
  data.frame(matrix(as.numeric(unlist(strsplit(sub("DP4=", replacement = "",x), ","))), 
                    nrow=length(strsplit(sub("DP4=", replacement = "",x), ",")), byrow = T))
}

SNP_reads <- function(x,y) { #Sum read number per allele ref and alt, add genotype status and put it in a dataframe
  data.frame(chr=x$chr, pos=x$pos, ref=x$ref, alt=x$alt, read_ref=rowSums(y[,1:2]), read_alt=rowSums(y[,3:4]),
             genotype=unlist(str_extract_all(unlist(strsplit(x$GT.ll,split=":")), "./.")))
}

inf10reads <- function(x) {x[(x$read_ref >= 3 & x$read_alt >= 3) & (x$read_ref + x$read_alt >= 10) & (x$genotype=="0/1" | x$genotype=="1/2" | x$genotype=="0/2"),] } 
# Remove SNPs covered by less than 10 reads & 1/1 genotype because they are homozygotes (diff from ref assembly)
#Keep only SNPs covered by atleast 3 reads and alt+ref = 10 reads

alt_freq <- function(x) {data.frame(x, alt_freq=1-(x$read_ref/rowSums(x[,5:6])))} # freq of alt allele

# SNPs chrX
chrX_DP4 <- extrac_DP4(filesSNPs[["chrX"]]$infos)
chrX_SNPreads <- extract_reads(chrX_DP4)
chrX_SNPs <- SNP_reads(filesSNPs[["chrX"]], chrX_SNPreads)
chrX_10rSNPs <- inf10reads(chrX_SNPs)
chrX_10rSNPs <- alt_freq(chrX_10rSNPs)

# SNPs chr7
chr7_DP4 <- extrac_DP4(filesSNPs[["chr7"]]$infos)
chr7_SNPreads <- extract_reads(chr7_DP4)
chr7_SNPs <- SNP_reads(filesSNPs[["chr7"]], chr7_SNPreads)
chr7_10rSNPs <- inf10reads(chr7_SNPs)
chr7_10rSNPs <- alt_freq(chr7_10rSNPs)

# SNPs chr8
chr8_DP4 <- extrac_DP4(filesSNPs[["chr8"]]$infos)
chr8_SNPreads <- extract_reads(chr8_DP4)
chr8_SNPs <- SNP_reads(filesSNPs[["chr8"]], chr8_SNPreads)
chr8_10rSNPs <- inf10reads(chr8_SNPs)
chr8_10rSNPs <- alt_freq(chr8_10rSNPs)

#Density plot indicating clonality regarding the Xi status
svg("clonality.svg", width = 17.5, height = 8) #saving in working directory
plot(density(chrX_10rSNPs$alt_freq),col ="blue", main ="SNPs expressions", xlab = "freq of alternative allele")
lines(density(chr7_10rSNPs$alt_freq), col = "pink")
lines(density(chr8_10rSNPs$alt_freq), col = "black")
legend("topright", legend=c("chrX", "chr7","chr8"), col = c("blue","pink","black"), lty = 1)
dev.off()

#Histogram of alt_SNPs frequencies confirmation of clonality
svg("hist_freq_alt_allele.svg", width = 17.5, height = 8)
par(mfrow=c(1,3))
hist(chrX_10rSNPs$alt_freq, main ="Allele freq expr SNPs chrX", xlab = "Alternative allele freq")
hist(chr7_10rSNPs$alt_freq, main ="Allele freq expr SNPs chr7", xlab = "Alternative allele freq")
hist(chr8_10rSNPs$alt_freq, main ="Allele freq expr SNPs chr8", xlab = "Alternative allele freq")
dev.off()

# SNPs densities

# load annotation data 
# if SNPs inferred from RNA-seq
#Concatenate chrX txn length
seqlevels(txdb) <- c("chrX") #Activates chrX only
TxnL_X <- transcriptLengths(txdb)
#Concatenate chr8 txn length
seqlevels(txdb) <- c("chr8") #Activates 8 only
TxnL_8 <- transcriptLengths(txdb)
#Concatenate chr7 txn length
seqlevels(txdb) <- c("chr7") #Activates 7 only
TxnL_7 <- transcriptLengths(txdb)

# if SNPs called from DNA-Seq
#txdb <- read.table(argmts[2], header = F, stringsAsFactors = F, sep = "\t", col.names = c("chr", "size"))
#Concatenate chrX txn length
#TxnL_X <- txdb[txdb$chr=="chrX",]
#Concatenate chr8 txn length
#TxnL_8 <- txdb[txdb$chr=="chr8",]
#Concatenate chr7 txn length
#TxnL_7 <- txdb[txdb$chr=="chr7",]


#chrX
dchrX <- data.frame(chr="chrX",SNPs=nrow(chrX_10rSNPs), chrSize=sum(TxnL_X$tx_len),
           SNPs_bp=nrow(chrX_10rSNPs)/sum(TxnL_X$tx_len), 
           SNPs_1kb=nrow(chrX_10rSNPs)/sum(TxnL_X$tx_len*1000))
#chr7
dchr7 <- data.frame(chr="chr7",SNPs=nrow(chr7_10rSNPs), chrSize=sum(TxnL_7$tx_len),
          SNPs_bp=nrow(chr7_10rSNPs)/sum(TxnL_7$tx_len), 
          SNPs_1kb=nrow(chr7_10rSNPs)/sum(TxnL_7$tx_len*1000))
#chr8
dchr8 <- data.frame(chr="chr8",SNPs=nrow(chr8_10rSNPs), chrSize=sum(TxnL_8$tx_len),
          SNPs_bp=nrow(chr8_10rSNPs)/sum(TxnL_8$tx_len), 
          SNPs_1kb=nrow(chr8_10rSNPs)/sum(TxnL_8$tx_len*1000))
          
dSNPs <- rbind(dchrX, dchr7, dchr8)

write.table(dSNPs, "SNPs_density.tsv", sep = "\t", dec=".", col.names = T, row.names = F)

# Generating bed files 
# SNPs 
infSNPs <- filesSNPs[["chrX"]][grep(pattern = ".*0/1|.*1/2", filesSNPs[["chrX"]]$GT.ll),] #Filtering non informative SNPs
X3mbSNP <- subset(infSNPs, infSNPs$pos > lbound & infSNPs$pos < rbound, select =c("chr", "pos", "ref","alt"))

# DpnII
dpn2sub <- subset(dpn2, dpn2$start > lbound & dpn2$end < rbound, select =c("chr","start","end"))

#bed files
dnp3xbed <- data.frame(chr=dpn2sub$chr, start=dpn2sub$start, end=dpn2sub$end, name="DpnII",score=0,strand=".", 
                       thickstart=dpn2sub$start, thickend=dpn2sub$end, itemRGB="51,153,255")# bleu
write.table(dnp3xbed, "XIC3mb_DpnII.bed", quote = F, sep = "\t", col.names = F, row.names = F, dec = ",")


X3mbSNPbed <- data.frame(chr=X3mbSNP$chr, start=X3mbSNP$pos, end=X3mbSNP$pos, name=X3mbSNP$alt,score=0,strand=".", 
                         thickstart=X3mbSNP$pos, thickend=X3mbSNP$pos, itemRGB="204,51,0")#rouge

write.table(X3mbSNPbed, "XIC3mb_SNPs.bed", quote = F, sep = "\t", col.names = F, row.names = F, dec = ",")
