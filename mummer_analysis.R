setwd("/home/emmanuel/Documents/Primates/Macaque/séquences/rheMac8/mummer/")

library("ggplot2")
library("svglite")
library("RColorBrewer")

show.coord <- read.table("XIC_rheMac8VScalJac3_mum.mcoords", sep = "\t", header = F)

colnames(show.coord) <- c("start_ref", "end_ref", "start_qry", "end_qry","length_alg_ref", "length_alg_qry",
                          "perc_id", "length_ref", "length_qry", "cov_ref","cov_qry","tag_ref",
                          "tag_qry")

# Replace with XIC sizes see lab book or below. 

#hg 38 1052756
#rhemac8 1035791
# panTro5 1120458
#calJac3 1220870

coord.qry <- 1035791-(1035791-show.coord$end_qry)

perc.id <- data.frame(coord.qry, perc_id=show.coord$perc_id, cov_qry = (1035791*(show.coord$cov_qry/100)))


ggplot(perc.id, aes(coord.qry, perc_id))+
  geom_line()+
  geom_smooth(method = "lm")+
  labs(x="rheMac8 XIC (bp)",y="% identity" )+
  theme_classic()
ggsave("XIC_rheMac8VScalJac3_percid.svg", width = 17.5, height = 8)

#Aligned block %coverage of query sequence
ggplot(perc.id, aes(coord.qry, cov_qry))+
  geom_line()+
  geom_hline(yintercept = mean(perc.id$cov_qry), color = "red", linetype="dashed", size=1.5)+
  labs(x="rheMac8 XIC (bp)",y="aligned blocks coverage (bp)" )+
  theme_classic()
ggsave("XIC_rheMac8VScalJac3_blockscov.svg", width = 17.5, height = 8)


