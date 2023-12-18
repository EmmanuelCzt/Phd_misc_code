#!/usr/bin/Rscript --vanilla
#############################
###     PARAMETERS	      ###
#############################
projectName="ncTx_recons"
organism="hg38" # "hg19" or "mm10"
expType="RNAseq" # "RNAseq"
annoType="stringtie" # Gencode or ncRNA


###################################
### DIRECTORIES and FILES       ###
###################################
workdir <- "~/Documents/Omics_dataset/RNASeq/nCoding_transcriptome/test"
res_folder_name <- file.path(workdir,"Stringtie_Gencode_without_overlapGenes");if(!file.exists(res_folder_name)){dir.create(res_folder_name)}

ProjectDir <- workdir;if(!file.exists(ProjectDir)){dir.create(ProjectDir)}
# resdir <- file.path(ProjectDir, res_folder_name);if(!file.exists(resdir)){dir.create(resdir)}

load("~/Desktop/Bladder2018/GTF/Gencode_v24_stringtie_2189_genes_17_7_17.RData")
load(file.path(res_folder_name, "exp.RData"))

###################################
###        Variables            ###
###################################

corThreshold <- 0.9
geneDistance <- 20000

###################################
###         Libraries           ###
###################################
library(pheatmap)
library(FactoMineR)
library(ggplot2)
library(factoextra)
library(scales)
library(devtools)
library("easyGgplot2")

#########################################
##  Find Groups of Genes to merge      ##
#########################################

#Create neighbour matrix
gt <- merged_all_gtf
rownames(gt) <- gt$gene_id


gtlinc <- gt
gtlinc <- gt[gt$gene_biotype %in% c("lincRNA", "lncRNA"),]; dim(gtlinc)
gtnonlinc <- gt[!(gt$gene_biotype %in% c("lincRNA", "lncRNA")),]; dim(gtnonlinc)

exp_linc <- exp[rownames(exp) %in% rownames(gtlinc),]; dim(exp_linc)
exp_linc <- exp_linc[order(match(rownames(exp_linc), rownames(gtlinc))),]; dim(exp_linc)

##############################################################
#   Histogram for neighbours genes correlation repartition   #
##############################################################

createCorrelationHistogram <- function(gt, exp, type, anno){
  gtlincStringtie <- gt
  gtlincStringtie <- gt[gt$gene_biotype %in% type,]; dim(gtlinc)
  exp_lincStringtie <- exp[rownames(exp) %in% rownames(gtlincStringtie),]; dim(exp_lincStringtie)
  exp_lincStringtie <- exp_lincStringtie[order(match(rownames(exp_lincStringtie), rownames(gtlincStringtie))),]; dim(exp_lincStringtie)

  cor_lincRNAsStringtie <- cor(t(exp_lincStringtie))
  #hist(cor_lincRNAsStringtie)

  cor_neighbours_lincs_Str <- NULL
  for (i in 1:(dim(cor_lincRNAsStringtie)[1] - 1)){
    cor_neighbours_lincs_Str <- c(cor_neighbours_lincs_Str, cor_lincRNAsStringtie[i,i+1])
  }
  corl <- data.frame(
    cor_nlS=cor_neighbours_lincs_Str[!(is.na(cor_neighbours_lincs_Str))]
  ) 
  ggplot(corl, aes(cor_nlS)) +
    geom_histogram(aes(y=..density.., fill = ..count..),
                   binwidth = 0.03, 
                   col="black") +
    labs(title=paste("Correlation scores in neighbours genes :", anno, "annotation")) +
    labs(x="Correlation score", y="Density") 
    ggsave(paste0(resdir,"/Histogram_correlation_repartition",anno,".pdf"), plot = last_plot(), width = 7, height = 8)
}



createCorrelationHistogram(gt,exp,c("lncRNA"), "StringTie")
createCorrelationHistogram(gt,exp,c("lincRNA"), "Gencode")
createCorrelationHistogram(gt,exp,c("lncRNA",  "lincRNA"), "Gencode + StringTie")



#Create a correlation matrix with all lincRNAs
#Genes are ordered by position
cor_lincRNAs <- cor(t(exp_linc))

cor_neighbours_lincs <- NULL
for (i in 1:(dim(cor_lincRNAs)[1] - 1)){
  cor_neighbours_lincs <- c(cor_neighbours_lincs, cor_lincRNAs[i,i+1])
}

hist(cor_lincRNAs,freq=FALSE)
qplot(cor_lincRNAs, geom="histogram")


chr1gt <- gtlinc[gtlinc$Chr == "chr1",]
cor_chr1 <- cor_lincRNAs[rownames(cor_lincRNAs) %in% rownames(chr1gt),colnames(cor_lincRNAs) %in% rownames(chr1gt)]
#Create heatmap with expresion correlation matrix
pdf(file.path(resdir, paste0("pheatmap_clusters_chr1.pdf")), onefile = FALSE)
pheatmap(cor_chr1, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

png(file.path(resdir, paste0("pheatmap_clusters_chr1.png")))
pheatmap(cor_chr1, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

#Correlation matrix in chromosome 18
chr18gt <- gtlinc[gtlinc$Chr == "chr18",]
cor_chr18 <- cor_lincRNAs[rownames(cor_lincRNAs) %in% rownames(chr18gt),colnames(cor_lincRNAs) %in% rownames(chr18gt)]
#Create heatmap with expresion correlation matrix
pdf(file.path(resdir, paste0("pheatmap_clusters_chr18.pdf")), onefile = FALSE)
pheatmap(cor_chr18, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

png(file.path(resdir, paste0("pheatmap_clusters_chr18.png")))
pheatmap(cor_chr18, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

#Correlation matrix in chromosome 17
chr17gt <- gtlinc[gtlinc$Chr == "chr17",]
cor_chr17 <- cor_lincRNAs[rownames(cor_lincRNAs) %in% rownames(chr17gt),colnames(cor_lincRNAs) %in% rownames(chr17gt)]
#Create heatmap with expresion correlation matrix
pdf(file.path(resdir, paste0("pheatmap_clusters_chr17.pdf")), onefile = FALSE)
pheatmap(cor_chr17, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

png(file.path(resdir, paste0("pheatmap_clusters_chr17.png")))
pheatmap(cor_chr17, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

#Function that create a list of couples of genes that have to be merged
findListCouplesToMerge <- function(rcor, corThreshold, geneDistance, rgt, listCouplesToMerge){
  for(i in 1:(dim(rcor)[1] - 1)){
    if(!is.na(rcor[i, i + 1])){
      if(rcor[i, i + 1] >= corThreshold & abs(rgt[i,]$End - rgt[i+1,]$Start) <= geneDistance){
        print(listCouplesToMerge)
        if(rgt[i,]$annotation != "Gencode_v24" | rgt[i+1,]$annotation != "Gencode_v24"){
          listCouplesToMerge <- rbind(listCouplesToMerge, cbind(rownames(rgt[i,]), rownames(rgt[i+1,])))
  
        }
      }
    } 
  }
  return(listCouplesToMerge)
}

gtPlus <- gtlinc[gtlinc$Strand == "+",]; dim(gtPlus)
cor_plus <- cor_lincRNAs[rownames(cor_lincRNAs) %in% rownames(gtPlus), colnames(cor_lincRNAs) %in% rownames(gtPlus)]
dim(cor_plus)
listCouplesToMerge <- findListCouplesToMerge(cor_plus, corThreshold, geneDistance, gtPlus, NULL)
dim(listCouplesToMerge)

gtMinus <- gtlinc[gtlinc$Strand == "-",]
cor_minus <- cor_lincRNAs[rownames(cor_lincRNAs) %in% rownames(gtMinus), colnames(cor_lincRNAs) %in% rownames(gtMinus)]
listCouplesToMerge <- findListCouplesToMerge(cor_minus, corThreshold, geneDistance, gtMinus, listCouplesToMerge)
dim(listCouplesToMerge)


genesToDelete <- listCouplesToMerge[listCouplesToMerge[,1] %in% listCouplesToMerge[,2],1]

listLongCouples <- cbind(listCouplesToMerge[!(listCouplesToMerge[,1] %in% genesToDelete),1], listCouplesToMerge[!(listCouplesToMerge[,2] %in% genesToDelete),2])
#listLongCouples[,3] <- NA

#Supprimer lignes avec gènes se trouvant dans la liste de genes a merger

merged_gt_linc <- gtlinc[!(rownames(gtlinc) %in% unique(c(listCouplesToMerge[,1], listCouplesToMerge[,2]))),]

dim(merged_gt_linc)


listCorrespondances <- NULL
for(j in 1:dim(listLongCouples)[1]){
  #lapply(listLongCouples, function(x){
  #print(x)
  x <- listCouplesToMerge[c(which(listCouplesToMerge[,1] == listLongCouples[j,1]):
                                                  which(listCouplesToMerge[,2] == listLongCouples[j,2])), 1]
  lengthcod <- sum(gtlinc$LengthCod[rownames(gtlinc) %in% x])
  listCluster <- c(x,listLongCouples[j,2])
  #dans listLongAll, dire si Stringtie ou non
  matchENSG <- grep("ENSG.", listCluster, value = TRUE)
  if(!(identical(matchENSG, character(0)))){
    #print(listCluster)
    newGene <- paste0("STRL", matchENSG[1])
    annot <- "merged_Stringtie_Gencode"
  } else{newGene <- paste0("STRL", j); annot <- "merged_Stringtie"}
  for(elmt in listCluster){
    bin <- data.frame(cbind(elmt, newGene))
    listCorrespondances <- rbind(listCorrespondances, bin)
  }
  #Create new lines with merged genes
  print(j)
  newLine <- gtlinc[rownames(gtlinc) == listLongCouples[j,1],]
  newLine$End <- gtlinc$End[rownames(gtlinc) == listLongCouples[j,2]]
  rownames(newLine) <- newGene
  newLine$gene_id <- newGene
  newLine$annotation <- annot
  newLine$LengthCod <- lengthcod
  merged_gt_linc <- rbind(merged_gt_linc, newLine)
}


merged_gt_linc <- merged_gt_linc[order(merged_gt_linc$Start),]
merged_gt_linc <- merged_gt_linc[order(merged_gt_linc$Chr),]

new_exp <- exp[rownames(exp) %in% rownames(merged_gt_linc),]; dim(new_exp)
new_exp <- new_exp[order(match(rownames(new_exp), rownames(merged_gt_linc))),]; dim(new_exp)

cor_newlincRNAs <- cor(t(new_exp))

#Merged data in chromosome 1
chr1newgt <- merged_gt_linc[merged_gt_linc$Chr == "chr1",]
cor_newchr1 <- cor_newlincRNAs[rownames(cor_newlincRNAs) %in% rownames(chr1newgt),colnames(cor_newlincRNAs) %in% rownames(chr1newgt)]

#Create heatmap with expresion correlation matrix
pdf(file.path(resdir, paste0("pheatmap_new_clusters_chr1_", corThreshold, "_", geneDistance, ".pdf")), onefile = FALSE)
pheatmap(cor_newchr1, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

png(file.path(resdir, paste0("pheatmap_new_clusters_chr1_", corThreshold, "_", geneDistance, ".png")))
pheatmap(cor_newchr1, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

#Merged data in chromosome 18
chr18newgt <- merged_gt_linc[merged_gt_linc$Chr == "chr18",]
cor_newchr18 <- cor_newlincRNAs[rownames(cor_newlincRNAs) %in% rownames(chr18newgt),colnames(cor_newlincRNAs) %in% rownames(chr18newgt)]

#Create heatmap with expresion correlation matrix
pdf(file.path(resdir, paste0("pheatmap_new_clusters_chr18_", corThreshold, "_", geneDistance, ".pdf")), onefile = FALSE)
pheatmap(cor_newchr18, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

png(file.path(resdir, paste0("pheatmap_new_clusters_chr18_", corThreshold, "_", geneDistance, ".png")))
pheatmap(cor_newchr18, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

#Merged data in chromosome 17
chr17newgt <- merged_gt_linc[merged_gt_linc$Chr == "chr17",]
cor_newchr17 <- cor_newlincRNAs[rownames(cor_newlincRNAs) %in% rownames(chr17newgt),colnames(cor_newlincRNAs) %in% rownames(chr17newgt)]

#Create heatmap with expresion correlation matrix
pdf(file.path(resdir, paste0("pheatmap_new_clusters_chr17_", corThreshold, "_", geneDistance, ".pdf")), onefile = FALSE)
pheatmap(cor_newchr17, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

png(file.path(resdir, paste0("pheatmap_new_clusters_chr17_", corThreshold, "_", geneDistance, ".png")))
pheatmap(cor_newchr17, cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, show_rownames = FALSE, legend = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()


#Replace name of genes merged by their new names in gff file
gff <- stringtie
gff2 <- gff$name
gff3 <- gff$name
dim(listCorrespondances)
for(i in 1:dim(listCorrespondances)[1]){
  print(i)
  print(listCorrespondances[i,1])
  print(listCorrespondances[i,2])
  #gff2 <- sapply(gff2, function(x) gsub(listCorrespondances[i,1],listCorrespondances[i,2], x))
  gff3 <- gsub(listCorrespondances[i,1],listCorrespondances[i,2], gff3)
}
gff$name <- gff3

#gff4 <- gff[gff$Type == "exon",]
#Reprendre les protein_coding pour faire le merging
#txdb <- makeTxDbFromGFF("~/Desktop/Bladder2018/Stringtie_Gencode_without_overlapGenes/gencode_v24_stringTie.gff",format="gtf")

merged_all_gtf <- rbind(merged_gt_linc, gtnonlinc)

save(gff, merged_all_gtf ,file=file.path(resdir,"mergedDatagff.RData"))
write.table(gff,file.path(resdir,"gencode_v24_stringTie.gff"),sep="\t",quote=F, row.names=FALSE)
write.table(merged_all_gtf, file=file.path(resdir,"gencode_v24_stringtie.gtf"),sep="\t", quote=F, row.names=FALSE)
write.table(listCorrespondances, file.path(resdir,"listCorrespondances.txt"),sep="\t",quote=F, row.names=FALSE)


#Create a pie chart
beforem <- data.frame(beforeafter= c("before merging", "before merging", "after merging", 
                                     "after merging", "after merging", "after merging"),
                                     annotation = c("Gencode_v24", "stringtie_specific", "Gencode_v24",
                                                    "stringtie_specific", "merged_Stringtie_Gencode",
                                                    "merged_Stringtie"), 
                      number = c(nrow(gtlinc[gtlinc$annotation == "Gencode_v24",]),
                                 nrow(gtlinc[gtlinc$annotation == "stringtie_overlap_MiTranscriptome",]) +
                                 nrow(gtlinc[gtlinc$annotation == "stringtie_specific",]),
                                 nrow(merged_gt_linc[merged_gt_linc$annotation == "Gencode_v24",]),
                                 nrow(merged_gt_linc[merged_gt_linc$annotation == "stringtie_overlap_MiTranscriptome",]) +
                                 nrow(merged_gt_linc[merged_gt_linc$annotation == "stringtie_specific",]),
                                 nrow(merged_gt_linc[merged_gt_linc$annotation == "merged_Stringtie_Gencode",]),
                                 nrow(merged_gt_linc[merged_gt_linc$annotation == "merged_Stringtie",])))
beforem$pos <- with(beforem, ave(number, annotation, FUN = function(x) cumsum(x) - 0.8*x))
library(dplyr)
beforem <- beforem %>% group_by(beforeafter) %>% mutate(pos=cumsum(number)-0.6*number)

#Create a pie chart that show annotation repartition before and after merging
ggplot(data = beforem) + 
  geom_bar(aes(x = "", y = number, fill = annotation), stat = "identity") +
  coord_polar(theta = "y") + theme(axis.title.x=element_blank(),
                                   axis.ticks.x=element_blank()) +
  facet_grid(beforeafter ~ ., scales = "fixed", as.table = FALSE) 
ggsave(paste0(resdir, "/before_after_merging_pie_charts.pdf"), plot = last_plot(), width = 7, height = 8)



