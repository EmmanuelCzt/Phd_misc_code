setwd("/home/emmanuel/Documents/Primates/marmoset/rmsk/")

library("ggplot2")
library("svglite")
library("RColorBrewer")

rmsk <- read.table("rmsk_calJac3.fa.out", sep = "\t")
colnames(rmsk) <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStarts", "genoEnd","genoLeft",
                    "strand","repName","repClass", "repFamily","repStart","repEnd","repLeft","id")
# write.table(rmsk, "", sep = "\t", col.names = T, quote = F)

dir.create("rmsk_calJac3")

rmsk.class <- as.vector(unique(rmsk$repClass))
rmsk.Family <- unique(rmsk$repFamily)

# Whole genome TE distribution class and family
perc.class <- NULL
perc.family <- NULL

for (i in 1:length(rmsk.class)) {
  perc.class[i] <- (nrow(rmsk[rmsk$repClass==rmsk.class[i],])/nrow(rmsk))*100
}

for (i in 1:length(rmsk.Family)) {
  perc.family[i] <- (nrow(rmsk[rmsk$repFamily==rmsk.Family[i],])/nrow(rmsk))*100
}

wg.class <- data.frame(rmsk.class, perc.class)
wg.class <- wg.class[order(wg.class$perc.class, decreasing = T),]
wg.family <- data.frame(rmsk.Family, perc.family)
wg.family <- wg.family[order(wg.family$perc.family, decreasing = T),]

# Whole genome class
ggplot(wg.class, aes(reorder(rmsk.class, -perc.class), perc.class))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=sprintf("%0.2f",perc.class)), vjust = -1)+ 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))+
  labs( x = "Class", y= "%")

ggsave("rmsk_class_calJac3.svg", path = "rmsk_calJac3/", width = 9, height = 7.5) #requires svglite

#Whole genome families
ggplot(wg.family, aes(reorder(rmsk.Family, -perc.family), perc.family))+
  geom_bar(stat = "identity", width = 0.7)+
  geom_text(aes(label=sprintf("%0.2f",perc.family)), vjust = 0, hjust = 0, angle = 45)+ 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))+
  labs( x = "Family", y= "%")

ggsave("rmsk_family_calJac3.svg", path = "rmsk_calJac3/", width = 11.6, height = 8) #requires svglite


# ChrX TE ditribution compared to the X chromosome 
rmsk.chrX <- subset(rmsk, genoName=="chrX", select = bin:id)

perc.class.chrX <- NULL
perc.family.chrX <- NULL

for (i in 1:length(rmsk.class)) {
  perc.class.chrX[i] <- ((nrow(rmsk.chrX[rmsk.chrX$repClass==rmsk.class[i],])/nrow(rmsk.chrX))*100)
}

i <- NULL

for (i in 1:length(rmsk.Family)) {
  perc.family.chrX[i] <- (nrow(rmsk.chrX[rmsk.chrX$repFamily==rmsk.Family[i],])/nrow(rmsk.chrX))*100
}

wg.class.chrX <- data.frame(rmsk.class, perc.class.chrX)
wg.family.chrX <- data.frame(rmsk.Family, perc.family.chrX)

ggplot(wg.class.chrX, aes(reorder(rmsk.class, -perc.class.chrX), perc.class.chrX))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=sprintf("%0.2f",perc.class.chrX)), vjust = -1)+ 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))+
  labs( x = "Class", y= "%")

ggsave("rmsk_class_chrX_calJac3.svg", path = "rmsk_calJac3/", width = 9, height = 7.5) #requires svglite


ggplot(wg.family.chrX, aes(reorder(rmsk.Family, -perc.family.chrX), perc.family.chrX))+
  geom_bar(stat = "identity", width = 0.7)+
  geom_text(aes(label=sprintf("%0.2f",perc.family.chrX)), vjust = 0, hjust = 0, angle = 45)+ 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))+
  labs( x = "Family", y= "%")

ggsave("rmsk_family_chrX_calJac3.svg", path = "rmsk_calJac3/", width = 11.6, height = 8)

# XIC
# See .md lab book for XIC coordinates in each species

rmsk.XIC <- read.table("rmsk_XIC_calJac3.fa.out", sep = "\t")
colnames(rmsk.XIC) <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStarts", "genoEnd","genoLeft",
                    "strand","repName","repClass", "repFamily","repStart","repEnd","repLeft","id")

rmsk.family.XIC <- unique(rmsk.XIC$repFamily)

perc.class.XIC <- NULL
perc.family.XIC <- NULL

for (i in 1:length(rmsk.class)) {
  perc.class.XIC[i] <- ((nrow(rmsk.XIC[rmsk.XIC$repClass==rmsk.class[i],])/nrow(rmsk.XIC))*100)
}

i <- NULL

for (i in 1:length(rmsk.family.XIC)) {
  perc.family.XIC[i] <- (nrow(rmsk.XIC[rmsk.XIC$repFamily==rmsk.family.XIC[i],])/nrow(rmsk.XIC))*100
}

wg.class.XIC <- data.frame(rmsk.class, perc.class.XIC)
wg.family.XIC <- data.frame(rmsk.family.XIC, perc.family.XIC)

ggplot(wg.class.XIC, aes(reorder(rmsk.class, -perc.class.XIC), perc.class.XIC))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=sprintf("%0.2f",perc.class.XIC)), vjust = -1)+ 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))+
  labs( x = "Class", y= "%")

ggsave("rmsk_class_XIC_calJac3.svg", path = "rmsk_calJac3/", width = 9, height = 7.5) #requires svglite


ggplot(wg.family.XIC, aes(reorder(rmsk.family.XIC, -perc.family.XIC), perc.family.XIC))+
  geom_bar(stat = "identity", width = 0.7)+
  geom_text(aes(label=sprintf("%0.2f",perc.family.XIC)), vjust = 0, hjust = 0, angle = 45)+ 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1))+
  labs( x = "Family", y= "%")

ggsave("rmsk_family_XIC_calJac3.svg", path = "rmsk_calJac3/", width = 11.6, height = 8)

#Density of RE along the XIC

RE.start <- data.frame((1220870-(66400555-rmsk.XIC$genoStarts)), rmsk.XIC$repFamily)
colnames(RE.start) <- c("RE_5", "repFamily") # See .md lab book for XIC sizes

#hg 38 1052756
#rhemac8 1035791
# calJac3 1120458
#calJac3 1220870

ggplot(RE.start, aes(x = RE_5))+
  geom_area(stat = "bin", binwidth = 5000)+
  scale_x_continuous(breaks=seq(0,1220870, 50000))+
  labs(x="XIC (bp)", y= "RE/5000bp")+
  theme_classic()
ggsave("rmsk_REdensity_XIC_calJac3.svg", path = "rmsk_calJac3/", width = 17.5, height = 8)

sup2 <- subset(RE.start, 
               repFamily==wg.family.XIC$rmsk.family.XIC[wg.family.XIC$perc.family.XIC > 5], 
               select = RE_5:repFamily) # subset families > 1.5%

ggplot(sup2, aes(x = RE_5, fill = repFamily))+
  geom_area(stat = "bin", binwidth = 30000)+
  scale_x_continuous(breaks=seq(0,1220870, 50000))+
  labs(x="XIC (bp)", y= "RE/30000bp")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()
ggsave("rmsk_REdensitysup5_XIC_calJac3.svg", path = "rmsk_calJac3/", width = 17.5, height = 8)
