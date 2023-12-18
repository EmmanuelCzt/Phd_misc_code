# Ct normalization

ct.long <- melt(ct)
colnames(ct.long)[4] <- "Primers"
ct.ae <- merge(ct.long,ae, by="Primers")
d.ct <- data.frame(ct.ae[,c(1:4)],ct=ct.ae$value,d=(ct.ae$AE^-ct.ae$value))
d.ct <- dcast(d.ct, Name ~ Primers)

do.norm <- function(x){x/d.ct$rhRPLPO}

dd.ct <- apply(X = d.ct[,-c(1,length(d.ct))],2,do.norm)
dd.ct <- data.frame(Name=d.ct$Name, dd.ct)

expr <- merge(ct[,c(1:3)], dd.ct, by="Name")

#Remove E10 because out
expr <- expr[-8,]


# Pooled gene expression

genot <- expr %>% group_split(Genot)

avrg <- lapply(genot, function(x){
  data.frame(Genot=unique(x$Genot), expr=apply(X = x[,-c(1:3)],2,mean),variable=colnames(x[,-c(1:3)]), row.names = NULL)
})

avrg <- lapply(avrg, function(x){data.frame(x,variable=rownames(x))})
avrg.variable <- do.call("rbind",avrg) %>% group_split(variable)
