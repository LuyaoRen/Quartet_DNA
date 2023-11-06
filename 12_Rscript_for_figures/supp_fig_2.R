
# load libraries
library(ggplot2)

###load data
filtered <- read.table('Quartet_DNA_BGI_SEQ2000_BGI_LCL5_1_20180518.filtered.vcf.info')
filtered <- filtered[,c(5,6,7,8)]
filtered$tag <- 'Filtered'

mendelian <- read.table('Quartet_DNA_BGI_SEQ2000_BGI_LCL5_1_20180518.mendelian.vcf.info')
mendelian <- mendelian[,c(5,6,7,8)]
mendelian$tag <- 'Mendelian Violation'

reference <- read.table('Quartet_DNA_BGI_SEQ2000_BGI_LCL5_1_20180518.ref.vcf.info')
reference <- reference[,c(5,6,7,8)]
reference$tag <- 'Reference'  

dat <- data.frame(rbind(filtered,mendelian,reference))
colnames(dat) <- c('DP','AF','GQ','MQ','Type')
dat$DP[which(dat$DP > 100)] <- 101
dat <- dat[which(dat$DP != 0),]
dat$AF <- as.numeric(dat$AF)

p <- ggplot(dat, aes(x=DP,color=Type, fill=Type)) + 
  geom_density(alpha=0.4)+
  #  scale_color_manual(values=c("#00bfc4","#f8766d"))+
  #  scale_fill_manual(values=c("#00bfc4","#f8766d"))+
  theme_light()+
  ylab('Density') + 
  xlab('Depth')+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

pdf('depth.pdf',height = 2.5,width = 4.5)
plot(p)
dev.off()


p <- ggplot(dat, aes(x=AF,color=Type, fill=Type)) + 
  geom_density(alpha=0.4)+
  #  scale_color_manual(values=c("#00bfc4","#f8766d"))+
  #  scale_fill_manual(values=c("#00bfc4","#f8766d"))+
  theme_light()+
  ylab('Density') + 
  xlab('Allele Frequency')+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

pdf('AF.pdf',height = 2.5,width = 4.5)
plot(p)
dev.off()


p <- ggplot(dat, aes(x=Type, y=GQ)) + 
  geom_violin(trim=FALSE, fill="gray")+
  labs(x="", y = "Genotype Quality")+
  stat_summary(fun=median, geom="point", size=2, color="red")+
  theme_classic()+
  theme(axis.text.x = element_text( color="black", 
                                    size=11, angle=45,hjust=1,vjust=1),
        axis.text.y = element_text( color="black", 
                                    size=11))

pdf('GQ.pdf',height = 3,width = 2.5)
plot(p)
dev.off()

p <- ggplot(dat, aes(x=Type, y=MQ)) + 
  geom_violin(trim=FALSE, fill="gray")+
  labs(x="", y = "Mapping Quality")+
  stat_summary(fun=median, geom="point", size=2, color="red")+
  theme_classic()+
  theme(axis.text.x = element_text( color="black", 
                                    size=11, angle=45,hjust=1,vjust=1),
        axis.text.y = element_text( color="black", 
                                    size=11))

pdf('MQ.pdf',height = 3,width = 2.5)
plot(p)
dev.off()

p <- ggplot(dat, aes(x=MQ,color=Type, fill=Type)) + 
  geom_density(alpha=0.4,position = "fill")+
  #  scale_color_manual(values=c("#00bfc4","#f8766d"))+
  #  scale_fill_manual(values=c("#00bfc4","#f8766d"))+
  theme_light()+
  ylab('Density') + 
  xlab('Mapping')+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))


####density plot

filtered$density <- densCols(filtered$V7, filtered$V8, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
mendelian$density <- densCols(mendelian$V7, mendelian$V8, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
reference$density <- densCols(reference$V7, reference$V8, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))

p <- ggplot(mendelian) +
  geom_point(aes(V7, V8, col = density), size = 0.5) +
  scale_color_identity() +
  ylab('Mapping Quality') + 
  xlab('Genotype Quality')+ 
  theme_bw()+theme(axis.text.x = element_text(size = 15),
                   axis.text.y = element_text(size = 15),
                   axis.title.x = element_text(size = 15),
                   axis.title.y = element_text(size = 15))
png('mendelian_mq_gq.png',res=120,height=400,width=500)
plot(p)
dev.off()

############################################################
# SV
############################################################

dat <- read.table('MNF_AF_pacbio_LCL5_minimap2_cutesv.txt',header=T)

dat$Var[which(dat$Var > 100)] <- 101

p <- ggplot(dat, aes(x=Var,color=class, fill=class)) + 
  geom_density(alpha=0.4)+
  #geom_histogram(aes(y=..density..), alpha=0.5, bins=50,
  #               position="identity")+
  #  scale_color_manual(values=c("#00bfc4","#f8766d"))+
  #  scale_fill_manual(values=c("#00bfc4","#f8766d"))+
  theme_light()+
  ylab('Density') + 
  xlab('Reads Supporting SV')+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

pdf('Reads_Supporting_SV.pdf',height = 2.5,width = 4.5)
plot(p)
dev.off()


p <- ggplot(dat, aes(x=AF,color=class, fill=class)) + 
  geom_density(alpha=0.4)+
  #geom_histogram(aes(y=..density..), alpha=0.5, bins=50,
  #               position="identity")+
  #  scale_color_manual(values=c("#00bfc4","#f8766d"))+
  #  scale_fill_manual(values=c("#00bfc4","#f8766d"))+
  theme_light()+
  ylab('Density') + 
  xlab('Allele Frequency')+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

pdf('SV_af.pdf',height = 2.5,width = 4.5)
plot(p)
dev.off()


dat <- read.table('MNF_MAPQ_pacbio_LCL5_minimap2_cutesv.txt',header=T)
p <- ggplot(dat, aes(x=class,y=V1, fill=class)) + 
  geom_boxplot(alpha=0.4)+
  #geom_violin(trim=FALSE)+
  #geom_histogram(aes(y=..density..), alpha=0.5, bins=50,
  #               position="identity")+
  #  scale_color_manual(values=c("#00bfc4","#f8766d"))+
  #  scale_fill_manual(values=c("#00bfc4","#f8766d"))+
  theme_light()+
  ylab('Mapping Quality') + 
  xlab('')+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))+
  theme(axis.text.x = element_text( color="black", 
                                    size=11, angle=45,hjust=1,vjust=1),
        axis.text.y = element_text( color="black", 
                                    size=11)) 

p <- ggplot(dat, aes(x=V1,color=class, fill=class)) + 
  geom_density(alpha=0.4)+
  #  scale_color_manual(values=c("#00bfc4","#f8766d"))+
  #  scale_fill_manual(values=c("#00bfc4","#f8766d"))+
  theme_light()+
  ylab('Density') + 
  xlab('Mapping Quality')+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

pdf('SV_MQ.pdf',height = 2.5,width = 4.5)
plot(p)
dev.off()


## density

filtered <- dat[dat$class == "Filter",]
mendelian <- dat[dat$class == "Mendelian_Violation",]
reference <- dat[dat$class == "Reference",]

filtered$density <- densCols(filtered$Var, filtered$AF, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
mendelian$density <- densCols(mendelian$Var, mendelian$AF, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
reference$density <- densCols(reference$Var, reference$AF, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))

p <- ggplot(reference) +
  geom_point(aes(Var, AF, col = density), size = 0.5) +
  scale_color_identity() +
  ylab('AF') + 
  xlab('Var')+ 
  xlim(c(0,100))+
  theme_bw()+theme(axis.text.x = element_text(size = 15),
                   axis.text.y = element_text(size = 15),
                   axis.title.x = element_text(size = 15),
                   axis.title.y = element_text(size = 15))
png('mendelian_mq_gq.png',res=120,height=400,width=500)
plot(p)
dev.off()
