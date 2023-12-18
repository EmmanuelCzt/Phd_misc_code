library(ggplot2)
library(dplyr)
library(stringr)
library(basicTrendline)

setwd("/home/emmanuel/Documents/qPCR/20210625_rhXICSetUp/")
raw <- read.table("/home/emmanuel/Documents/qPCR/20210625_rhXICSetUp/2021-06-25 rhXICSetUp.csv", header = T, sep = "\t", dec = ",", stringsAsFactors = F)

#SD >0.5

sdhigh <- subset(raw, Ct.SD>=0.5)
write.table(sdhigh,"SDhigh.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
raw[raw$Ct.SD>=0.5,"warning"] <- "yes"

#More than one PCR product
PCR2 <- subset(raw,!(Tm2=="" & Tm3==""))
write.table(PCR2, "NonspePCRproduct.tsv",sep = "\t", quote = F, col.names = T, row.names = F)
raw[!(raw$Tm2=="" & raw$Tm3==""),"warning"] <- "yes"
raw[(raw$Tm2=="" & raw$Tm3==""),"warning"] <- "no"
#Remove RTneg
raw <- raw[raw$Sample.Name!="RTNEG",]

#Sample name
target <- unique(raw$Target.Name)

paramlm <- list()

dir.create("fig")
dflm <- data.frame(matrix(nrow=length(target),ncol=8))
colnames(dflm) <- c("a","b","R","AE","eff","Tm1","Tm1sd","Name")

for(i in 1:length(target)){
  paramlm <- lm(formula = as.numeric(CT)~log(as.numeric(Sample.Name)), data = subset(raw, Target.Name==target[i]))
  
  dflm[i,1] <- coef(paramlm)[2]
  dflm[i,2] <- coef(paramlm)[1]
  dflm[i,3] <- summary(paramlm)$r.squared
  dflm[i,4] <- 10^(-(1/(coef(paramlm)[2]*log(10))))
  dflm[i,5] <- 10^(-(1/(coef(paramlm)[2]*log(10))))/2
  dflm[i,8] <- target[i]
  
  if(raw$warning[i]=="yes"){
    ggplot(subset(raw, Target.Name==target[i]), aes(x=as.numeric(Sample.Name),y=as.numeric(CT)))+
      geom_point(colour="red")+
      geom_smooth( method="lm")+
      ggtitle(paste(target[i]," - ","y=",round(dflm$a[i],3)," x ln(x) + ",round(dflm$b[i],3)," - R2=",round(dflm$R[i],3), " - AE=",round(dflm$AE[i],3),sep = ""))+
      scale_y_continuous(limits = c(20.0,30.0))+
      scale_x_log10()+
      ylab("CT")+
      xlab("log(cDNA dilution)")
    ggsave(paste("fig/",target[i],".jpeg", sep = ""),device = "jpeg")
  }else {
  ggplot(subset(raw, Target.Name==target[i]), aes(x=as.numeric(Sample.Name),y=as.numeric(CT)))+
    geom_point()+
    geom_smooth( method="lm")+
    ggtitle(paste(target[i]," - ","y=",round(dflm$a[i],3)," x ln(x) + ",round(dflm$b[i],3)," - R2=",round(dflm$R[i],3), " - AE=",round(dflm$AE[i],3),sep = ""))+
    scale_y_continuous(limits = c(15.0,35.0))+
    scale_x_log10()+
    ylab("CT")+
    xlab("log(cDNA dilution)")
  ggsave(paste("fig/",target[i],".jpeg", sep = ""),device = "jpeg")
  }
}

for(i in 1:length(target)){
  dflm[i,"Tm1"] <- mean(as.numeric(subset(raw, Target.Name==target[i])$Tm1))
  dflm[i,"Tm1sd"] <- sd(as.numeric(subset(raw, Target.Name==target[i])$Tm1))
}

write.table(dflm,"qPCR_primers_stats.tsv",sep = "\t", dec = ".", col.names = T, quote = F, row.names = F)
