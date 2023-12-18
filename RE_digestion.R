### Load the needed libraries
# Biostrings is needed for pattern identification
# BSgenome.Mmusculus.UCSC.mm10 is the mouse genome
# plyr and reshape2 are needed for manipulating the data format
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.hg38)
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)

# Identify the MspI recongnition sites for each chromosomal entry
# Generate a dataframe with the length of MspI digested fragments
mdf=data.frame();
for (i in seq_along(Mmusculus)){
  #print(paste("Processing ",seqnames(Mmusculus)[i], sep=""))
  m<-matchPattern("CCGG", Mmusculus[[i]])
  starts<-start(gaps(m))
  ends<-end(gaps(m))
  temp_df<-data.frame(start=starts-4,end=ends,chr=seqnames(Mmusculus)[i]) #actually end = ends
  temp_df$start<-replace(temp_df$start, temp_df$start == -3, 0)
  temp_df<-temp_df[c("chr","start","end")]
  mdf<-rbind(mdf,temp_df)
}

# Extract the digested fragment length
# Keep only the fragments in the range of 40-600 bp
mdf$width=mdf$end-mdf$start
ml<-mdf[mdf$width>39&mdf$width<601,]
counts<-ddply(ml,.(width), nrow)
# Create a plot of the frequency of the fragment lengths (y-axis is logarithmic)
p<-ggplot(counts,aes(x=width, y=V1))+geom_line()
p+scale_y_continuous(trans=log2_trans())