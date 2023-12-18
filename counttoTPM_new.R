library(GenomicFeatures)
countsToTPM <- function(raw.counts, tx.model,outfile){
  #Extract tx model
  gtf <- makeTxDbFromGFF(tx.model,format = "gtf")
  #Extract Tx model per exon
  exonlist <- exonsBy(gtf, by="gene")
  #Get tx sizes
  exonic.gene.sizes <- lapply(exonlist,function(x){sum(width(reduce(x)))})
  transcripts.length <- data.frame(row.names = names(exonic.gene.sizes), tx_length=matrix(unlist(exonic.gene.sizes)))
  #Merge with raw counts
  counts.tx.sizes <- merge(raw.counts, transcripts.length, by="row.names")
  # compute RPK  
  RPK <- sapply(counts.tx.sizes[,2:(length(counts.tx.sizes)-1)], function(x){x*1000/counts.tx.sizes$tx_length})
  #compute TPM
  TPM <- sapply(as.data.frame(RPK), function(x){(x/sum(x))*1000000})
  rownames(TPM) <- counts.tx.sizes$Row.names
  write.table(TPM, file = paste(outfile,".TPM", sep = ""), col.names = T, row.names = T, sep = "\t")
  return(TPM)
}