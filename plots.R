##RT-qPCR Barplot

ggplot()+
  geom_bar(data=avrg, aes(x = Genot, y=avrg[,genenames[1]], fill=Genot),
           stat = "identity",width = .5, 
           position=position_dodge(width = 0.5), colour="black", 
           size=1.5, fill=c((pal_material("blue-grey")(10))[c(10,6,1)]))+
  geom_point(data=expr, aes(x=Genot,y=expr[,genenames[1]]), size=3, colour="black")+
  ylab("Relative expression to RPLP0 mRNA levels")+
  xlab("")+
  scale_y_continuous( expand = expansion(mult = c(0,0.1)))+
  ggtitle(genenames[1])+
  theme_pubr()+
  theme(axis.text  = element_text(size=15,face = "bold", colour = "black", family = "Arial"),
        axis.text.x = element_text(angle = .45, family = "Arial"),
        axis.ticks =element_line(size=1, colour = "black"), axis.ticks.length = unit(.2,"cm"),
        axis.line = element_line(size=1, colour = "black"),
        axis.title = element_text(size=15, family = "Arial"),
        aspect.ratio = 2/1,
        legend.position = "none")


##percent stacked barplot RNA-FISH
ggplot()+
  geom_bar(data=xist,aes(x = Genot,y=counts, fill=pattern),position="fill",stat = "identity", colour="black", width =.6,size=1)+
  geom_text(data = del.lab, aes(x = x, y=y, label=labs),size=6)+
  scale_fill_manual(values=pal_material("grey")(10)[c(10,5,7)])+
  scale_y_continuous( expand = expansion(mult = c(0,0.1)),breaks = seq(0,1,0.2))+
  xlab("")+
  ylab("Observed cell fraction")+
  theme_pubr()+
  theme(axis.text  = element_text(size=15,face = "bold", colour = "black", family = "Arial"),
        axis.ticks =element_line(size=1, colour = "black"), axis.ticks.length = unit(.2,"cm"),
        axis.line = element_line(size=1, colour = "black"),
        axis.title = element_text(size=15, family = "Arial"),
        strip.text = element_text(size=15,face = "bold", colour = "black", family = "Arial"),
        legend.position = "right")