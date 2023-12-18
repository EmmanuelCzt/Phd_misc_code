#!/usr/bin/env Rscript

#argmts <- commandArgs(trailingOnly = T) # count file = arg 1 & gtf file = arg 2 

library(GenomicFeatures)

raw.counts <- read.table(snakemake@input[["count"]], sep = "\t", header = F, dec = ".")
raw.counts <- raw.counts[1:(dim(raw.counts)[1]-5),] # remove 5 last lines that are not usefull 
colnames(raw.counts) <- c("gene_id","counts")

# Import gene model
genemodel <- makeTxDbFromGFF(snakemake@input[["gtf"]], format = "gtf")

# Extract transcript length
exonslist <- exonsBy(genemodel, by="gene") #extract exons list by genes
exonic.gene.sizes <- lapply(exonslist,function(x){sum(width(reduce(x)))}) #create a list of genes lengths


transcripts.length <- data.frame(gene_id=names(exonic.gene.sizes), tx_length=matrix(unlist(exonic.gene.sizes)))

# Merge transcripts info and counts
rawcounts.to.TPM <- merge(transcripts.length, raw.counts, by = "gene_id")

# Formula TPM from wagner et al 2012
#Read per kilobase
rawcounts.to.TPM[,4] <- rawcounts.to.TPM$counts/rawcounts.to.TPM$tx_len
colnames(rawcounts.to.TPM)[4] <- c("RPK")

#per million scalling factor
rawcounts.to.TPM[,5] <-(rawcounts.to.TPM$RPK/(sum(rawcounts.to.TPM$RPK)))*1000000
colnames(rawcounts.to.TPM)[5] <- c("TPM")

#write.table(raw.counts, file=paste(snakemake@output[["count"]], ".raw", sep=""), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(rawcounts.to.TPM, file=snakemake@output[["out"]], sep = "\t", col.names = T, row.names = F, quote = F)