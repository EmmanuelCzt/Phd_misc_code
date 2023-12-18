########################################
#                                      #
#   Normalize HTSeq using TPM method   #
#          Louis Chauviere             #
#            22/08/2018                #
#                                      #
########################################

library(dplyr)
setwd("~/Bureau/single_cell/giustacchini/htseq2/")
listFiles <- list.files()
tissue_type <- "HSC"

#Create a HTSeq count table with all samples
counts.df <- NULL
for (i in c(1:length(listFiles))){
  countsFile <- read.table(listFiles[i])
  colnames(countsFile) <- c("gene_ID", unlist(strsplit(listFiles[i], "_"))[1])
  counts.df <- cbind(counts.df, countsFile[,2])
}
rownames(counts.df) <- countsFile[,1]
colnames(counts.df) <- matrix(unlist(strsplit(listFiles, ".", fixed = TRUE)), ncol=2, byrow = TRUE)[,1] #Find colnames with files names
counts.df <- counts.df[1:(dim(counts.df)[1]-5),] #delete last rows that are not genes
counts.df <- data.frame(counts.df)
counts.df <- cbind(rownames(counts.df), counts.df)
colnames(counts.df)[1] <- "Name"

#load bed file
gtf <- read.csv(file = "~/Documents/annotation/without_header_gencode.v29.annotation_XACT_T113_mitranscriptome.gff", sep = "\t", header = FALSE)

#calculate gene length using the gtf file and only the exons
workingTable <- gtf[gtf$V3 == "exon",]
library("GenomicFeatures")
txdb <- makeTxDbFromGFF(file = "~/Documents/annotation/gencode.v29.annotation_XACT_T113_mitranscriptome.gff", format = "gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene") #create a list of exons per genes
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))}) #create a list of genes lengths

gene.length <- as.data.frame(exonic.gene.sizes)
g <- matrix(gene.length)
g <- as.data.frame(g)
gene.length <- cbind(colnames(gene.length), g)


colnames(gene.length) <- c("Name", "gene_length")

library(dplyr)
counts.df <- inner_join(counts.df, gene.length, by="Name")


#####################################
#                                   #
#   Function for TPM calculation    #
#                                   #
#####################################

counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  return(tpm)
}

#Create a list of gene lengths 
featureLength <- unlist(counts.df$gene_length)
meanFragmentLength <- rep.int(51, (dim(counts.df)[2] - 2)) #mean fragment length is 350bp -> cf GSE69239


cdf <- sapply(counts.df[,2:(length(counts.df) -1)], as.numeric)
counts <- as.matrix(cdf)
rownames(counts) <- counts.df$Name

#Counts TPM from HTSeq counts
TPM.from.HTSeq <- counts_to_tpm(counts, featureLength, meanFragmentLength)
trans.TPM.from.HTSeq <- t(TPM.from.HTSeq)
TPM.from.HTSeq <- data.frame(TPM.from.HTSeq)


#boxplot(TPM.from.HTSeq, ylim=c(0,10)) #Normalization seems to be ok


#library(FactoMineR)
#PCA(X = TPM.from.HTSeq)




#load metadata
metadata <- read.csv("../Giustacchini_metadata.txt", 
                     header=T, sep = "\t")
#metadata$cell_type <- c(rep(tissue_type, dim(metadata)[1]))
#tumor_stage = cell_line

colnames(gtf) <- c("chr", "source", "type", "start", "end", "w", "strand", "y", "name")

#annotate count table
split<-function(x) {unlist(strsplit(x,split=";"))}
split2<-function(x) {unlist(strsplit(x,split=" "))}

workingTable <- gtf[gtf$type == "gene",]
workingTable$chr <- as.character(workingTable$chr)
workingTable$strand <- as.character(workingTable$strand)
workingTable$details <- sapply(as.character(workingTable$name), function(x) split(x),simplify = TRUE)

workingTable$gene_name3<-sapply(workingTable$details, function(x) x[grep("gene_name",x)],simplify=TRUE)
workingTable$gene_name<-sapply(workingTable$gene_name3,function(x) split2(x)[3],simplify=TRUE)

workingTable$gene_id3<-sapply(workingTable$details, function(x) x[grep("gene_id",x)],simplify=TRUE)
workingTable$gene_id<-sapply(workingTable$gene_id3,function(x) split2(x)[2],simplify=TRUE)

workingTable$gene_biotype3<-sapply(workingTable$details, function(x) x[grep("gene_type",x)],simplify=TRUE) #gene_type/
workingTable$gene_biotype<-sapply(workingTable$gene_biotype3,function(x) split2(x)[3],simplify=TRUE)



bedAll <- cbind(workingTable$chr, workingTable$start, workingTable$end, workingTable$strand, 
                workingTable$gene_name, workingTable$gene_id, workingTable$gene_biotype)
colnames(bedAll) <- c("chr", "start", "end", "strand", "gene_name", "gene_id", "gene_biotype")
bedAll <- data.frame(bedAll)

write.table(TPM.from.HTSeq, "../Giustacchini_554_cells_TPM_hg38.txt",sep=";",quote=F, row.names=TRUE)
save(TPM.from.HTSeq, file = "../Giustacchini_554_cells_TPM_hg38.RData")

type_TPM <- TPM.from.HTSeq

metadata <- metadata[metadata$Run %in% colnames(type_TPM),]

#colnames(type_TPM) <- metadata$tumor_stage

#annotate TPM count file
ann.TPM.from.HTSeq <- cbind(rownames(type_TPM), type_TPM)
ann.TPM.from.HTSeq$gene_id <- ann.TPM.from.HTSeq$`rownames(type_TPM)`
ann.TPM.from.HTSeq <- ann.TPM.from.HTSeq[,2:dim(ann.TPM.from.HTSeq)[2]]
ann.TPM.from.HTSeq$gene_id <- as.factor(ann.TPM.from.HTSeq$gene_id)

ann.TPM.from.HTSeq <- inner_join(ann.TPM.from.HTSeq, bedAll, by="gene_id")
dim(ann.TPM.from.HTSeq)
#rownames(ann.TPM.from.HTSeq) <- ann.TPM.from.HTSeq$gene_name


ann.TPM.from.HTSeq <- t(ann.TPM.from.HTSeq)
ann.TPM.from.HTSeq <- data.frame(ann.TPM.from.HTSeq)
#colnames(ann.TPM.from.HTSeq) <- ann.TPM.from.HTSeq[560,]
#ann.TPM.from.HTSeq <- ann.TPM.from.HTSeq[order(rownames(ann.TPM.from.HTSeq)),]

#barplot for XIST for each sample
png("../barplot_genes.png")
par(mfrow=c(3,2))

xist <- ann.TPM.from.HTSeq[ann.TPM.from.HTSeq$gene_name == "XIST",1:dim(TPM.from.HTSeq)[2]]
barplot(as.matrix(xist), las=2, cex.names = 0.8, main = "XIST", axisnames = FALSE)

ftx <- ann.TPM.from.HTSeq[ann.TPM.from.HTSeq$gene_name == "FTX",1:dim(TPM.from.HTSeq)[2]]
barplot(as.matrix(ftx), las=2, cex.names = 0.8, main = "FTX", axisnames = FALSE)

jpx <- ann.TPM.from.HTSeq[ann.TPM.from.HTSeq$gene_name == "JPX",1:dim(TPM.from.HTSeq)[2]]
barplot(as.matrix(jpx), las=2, cex.names = 0.8, main = "JPX", axisnames = FALSE)

gapdh <- ann.TPM.from.HTSeq[ann.TPM.from.HTSeq$gene_name == "GAPDH",1:dim(TPM.from.HTSeq)[2]]
barplot(as.matrix(gapdh), las=2, cex.names = 0.8, main = "GAPDH", axisnames = FALSE)

xact <- ann.TPM.from.HTSeq[ann.TPM.from.HTSeq$gene_name == "XACT",1:dim(TPM.from.HTSeq)[2]]
barplot(as.matrix(xact), las=1, cex.names = 0.8, main = "XACT")

t113 <- ann.TPM.from.HTSeq[ann.TPM.from.HTSeq$gene_name == "T113",1:dim(TPM.from.HTSeq)[2]]
barplot(as.matrix(t113), las=1, cex.names = 0.8, main = "T113")

dev.off()

ann_counts.df <- inner_join(counts.df, bedAll, by=c("Name"="gene_id"))
bedY<- ann_counts.df[ann_counts.df$chr == "chrY",2:(dim(TPM.from.HTSeq)[2] + 1)]
#bedY[] <- lapply(bedY, function(x)
#  as.numeric(levels(x))[x])

boxplot(bedY)#, ylim = c(0,1))
colSums(bedY)
