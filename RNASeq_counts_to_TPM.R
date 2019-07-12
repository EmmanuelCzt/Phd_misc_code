#!/usr/bin/env Rscript

argmts <- commandArgs(trailingOnly = T) # count file = arg 1 & gtf file = arg 2 

library("GenomicFeatures")

raw.counts <- read.table(argmts[1], sep = "\t", header = F, dec = ".")
raw.counts <- raw.counts[1:(dim(raw.counts)[1]-5),] # remove 5 last lines that are not usefull 
colnames(raw.counts) <- c("gene_id","counts")

# Import gene model
genemodel <- makeTxDbFromGFF(argmts[2], format = "gtf")

# Extract transcript length
transcripts.info <- transcriptLengths(genemodel)
transcripts.length <- transcripts.info[,c(3,5)] # Keep gene_id and length

# Merge transcripts info and counts
rawcounts.to.TPM <- merge(transcripts.length, raw.counts, by = "gene_id")

# Formula TPM from wagner et al 2012
#Read per kilobase
rawcounts.to.TPM[,4] <- rawcounts.to.TPM$counts/rawcounts.to.TPM$tx_len
colnames(rawcounts.to.TPM)[4] <- c("RPK")

#per million scalling factor
rawcounts.to.TPM[,5] <-(rawcounts.to.TPM$RPK/(sum(rawcounts.to.TPM$RPK)))*1000000
colnames(rawcounts.to.TPM)[5] <- c("TPM")

write.table(raw.counts, argmts[1], sep = "\t", col.names = T, row.names = F, quote = F)
write.table(rawcounts.to.TPM, file=paste(argmts[1], ".TPM", sep="") , sep = "\t", col.names = T, row.names = F, quote = F)