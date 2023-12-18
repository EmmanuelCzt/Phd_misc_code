setwd("/home/emmanuel/Documents/capture_HiC/large_cellCult/try2_3cellFate_markers/")
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

data <- read.table("2020_08_07_hEScH2-3_rhEScH2.txt",header = T, sep = "\t", stringsAsFactors = F, dec = ",")
data <- data[,-6]
AE <- read.table("/home/emmanuel/Documents/Bio_Mol/Oligos/AE.csv", header = F, sep = ";", stringsAsFactors = F, dec = ","
                 , row.names = 1)
colnames(AE) <- "AE"

#filter sd >=.5
sdhigh <- data[data$Ct.SD >= .5,]
write.table(sdhigh,"sdabove05.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#data <- data[data$Ct.SD < .5,]

#check Tm
ggplot(data, aes(x = Sample.Name, y=as.numeric(Tm1)))+
  geom_point()+
  coord_flip()+
  ylab("Tm °C")+
  xlab("Samples")+
  facet_wrap(~Target.Name)

ggsave("Tm_check.pdf", device = "pdf", width = 10, height = 17, units = "cm")

#remove duplicates
data <- data[seq(from=1, to=nrow(data), by=2),]

#Build ct table
samples <- unique(data$Sample.Name)
ct <- data.frame(matrix(nrow = length(unique(data$Target.Name)), 
                        ncol = length(unique(data$Sample.Name))),
                 row.names = unique(data$Target.Name))


for(i in 1:length(samples)){
  colnames(ct)[i] <- samples[i]
  ct[,i] <- as.numeric(data[data$Sample.Name==samples[i], 4])
}

#add AE
ct <- merge(ct, AE, by="row.names")
ct <- data.frame(ct, row.names = 1)

jpeg("CtrhRPL13A.jpg")
par(mfrow=c(1,2), mar=c(8,5,5,5))
barplot(as.matrix(ct["rhRPL13A", c(1:length(ct)-1)]),las =3, ylab = "Ct", main = "rhRPL13A Ct", las=2)
barplot(as.matrix(ct["hbACT", c(1:length(ct)-1)]),las =3, ylab = "Ct", main = "hbACT Ct", las=2)
dev.off()

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
#gene expr
gnexpr <- data.frame(matrix(nrow = length(unique(data$Target.Name)), 
                            ncol = length(unique(data$Sample.Name))),
                     row.names = rownames(ct))


for(i in 1:length(unique(data$Target.Name))){
  gnexpr[i,] <- ct[i,"AE"] ** -(ct[i,])
}
colnames(gnexpr) <- colnames(ct)[1:length(ct)-1]

#rhesus
gnexpr_rh <- gnexpr[grep("rh.*",x = row.names(gnexpr)),]
norm_rh <- data.frame(matrix(nrow = nrow(gnexpr_rh), 
                             ncol = length(unique(data$Sample.Name))),
                      row.names = row.names(gnexpr_rh))

for(i in 1:length(samples)){
  colnames(norm_rh)[i] <- samples[i]
  norm_rh[,i] <- gnexpr_rh[,i]/c(gnexpr_rh["rhRPL13A",i])
}
norm_rh <- norm_rh[, -c(1:4,11,12)]
#norm_rh[is.na(norm_rh)] <- 0
mn_rh <- data.frame(matrix(ncol = nrow(norm_rh), nrow = length(norm_rh)))
sd_rh <- data.frame(matrix(ncol = nrow(norm_rh), nrow = length(norm_rh)))

for(i in 1:nrow(norm_rh)){
  for(j in seq(1,length(norm_rh), by=2)){
    mn_rh[j,i] <- mean(as.numeric(norm_rh[i,j:c(j+1)]))
    sd_rh[j,i] <- sd(as.numeric(norm_rh[i,j:c(j+1)]), na.rm = TRUE)
  }
}

mn_rh <- mn_rh[-seq(0,length(norm_rh), by=2),]
sd_rh <- sd_rh[-seq(0,length(norm_rh), by=2),]
colnames(mn_rh) <- row.names(norm_rh)
colnames(sd_rh) <- row.names(norm_rh)

mn_rh$samples <- factor(c("rhES_B6","rhES_100","rhES_150"), 
                        levels=c("rhES_B6","rhES_100","rhES_150"))
sd_rh$samples <- factor(c("rhES_B6","rhES_100","rhES_150"),
                        levels=c("rhES_B6","rhES_100","rhES_150"))


rh <- list()
for(i in 1:c(length(mn_rh)-1)){
  rh[[i]] <- data.frame(mn=mn_rh[,i], sd=sd_rh[,i], smp=mn_rh$samples)
  names(rh)[i] <- colnames(mn_rh)[i]
}

for(i in seq_along(rh)){
  plot <- ggplot(rh[[i]], aes(x = smp, y=mn))+
          geom_bar(stat = "identity")+
          geom_errorbar( aes(x=smp, ymin=mn-sd, ymax=mn+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)+
          ggtitle(names(rh[i]))+
          xlab("Samples")+
          ylab("Expression relative to rhRPL13A mRNA")+
          theme_pubr()
  #print(plot)
  ggsave(plot, filename = paste("rhES_cHic_try2", names(rh[i]), ".jpg", sep=""), device = "jpg")
}

#human
gnexpr_hm <- gnexpr[grep("rh.*",x = row.names(gnexpr), invert = T),]
norm_hm <- data.frame(matrix(nrow = nrow(gnexpr_hm), 
                             ncol = length(unique(data$Sample.Name))),
                      row.names = row.names(gnexpr_hm))

for(i in 1:length(samples)){
  colnames(norm_hm)[i] <- samples[i]
  norm_hm[,i] <- as.numeric(gnexpr_hm[,i]/c(gnexpr_hm["hbACT",i]))
}
norm_hm <- data.frame(t(norm_hm[, c(1:4)]))

norm_hm$samples <- factor(rownames(norm_hm),levels = rownames(norm_hm))


hm <- list()
for(i in 1:c(length(norm_hm)-1)){
  hm[i] <- list(data.frame(expr=norm_hm[,i], smp=norm_hm$samples))
  names(hm)[i] <- colnames(norm_hm)[i]
}

for(i in seq_along(hm)){
  plot <- ggplot(hm[[i]], aes(x = smp, y=expr))+
    geom_bar(stat = "identity")+
    ggtitle(names(hm[i]))+
    xlab("Samples")+
    ylab("Expression relative to hbACT mRNA")+
    theme_pubr()
  #print(plot)
  ggsave(plot, filename = paste("hES_cHic_try2-3", names(hm[i]), ".jpg", sep=""), device = "jpg")
}
