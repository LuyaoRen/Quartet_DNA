library(reshape2)
library(ggplot2)
require("ggrepel")
library(cowplot)
dat <- read.table('nctr_precision_recall.txt',header=T)
dat$f1 <- 2*dat$precision*dat$recall/(dat$precision+dat$recall)

ref <- data_summary(dat, varname="f1", 
                    groupnames=c("mapper",'caller','site','gatk','type'))

dat2 <- read.table('nctr_quartet.txt',header=T)
quartet <- data_summary(dat2, varname="mendelian", 
                        groupnames=c("mapper",'caller','site','gatk','type'))
#scatter plot
nctr <- ref
nctr$mendelian <- quartet$mendelian
nctr$mendelian_sd <- quartet$sd
snv <- nctr[which(nctr$type == 'SNV'),]
p <- ggplot(snv, aes(x=f1, y=mendelian,shape=mapper, color=caller)) + geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))

pdf('nctr_snv_scatter.pdf',height=4,width=5)
plot(p)
dev.off()

indel <- nctr[which(nctr$type == 'INDEL'),]
p <- ggplot(indel, aes(x=f1, y=mendelian,shape=mapper, color=caller)) + geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))
pdf('nctr_indel_scatter.pdf',height=4,width=5)
plot(p)
dev.off()


snv <- dat[which(dat$type == 'INDEL'),]
p <- ggplot(snv, aes(x=precision, y=recall,shape=mapper, color=caller)) + geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))
pdf('indel_precision_recall.pdf',height=4,width=5.5)
plot(p)
dev.off()



 #barplot
ref <- data_summary(dat, varname="f1", 
                    groupnames=c('caller','type'))
snv_ref <- ref[grep('SNV',ref$type),]
p <- ggplot(snv, aes(x=caller, y=f1)) + 
  geom_bar(stat="identity", color='black',fill="#f16745", position=position_dodge()) +
  geom_errorbar(aes(ymin=f1-sd, ymax=f1+sd), width=.2,
                position=position_dodge(.9))+
  coord_cartesian(ylim = c(0.4, 1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))
pdf('nctr_indel_f1_barplot.pdf',height = 4,width=4)
plot(p)
dev.off()


quartet <- data_summary(dat2, varname="mendelian", 
                        groupnames=c('caller','type'))
snv_quartet <- quartet[grep('SNV',quartet$type),]
p <- ggplot(snv, aes(x=caller, y=mendelian)) + 
  geom_bar(stat="identity", color='black',fill="#177e89", position=position_dodge()) +
  geom_errorbar(aes(ymin=mendelian-sd, ymax=mendelian+sd), width=.2,
                position=position_dodge(.9))+
  coord_cartesian(ylim = c(0.9, 1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))
pdf('nctr_snv_mendelian_barplot.pdf',height = 4,width=4)
plot(p)
dev.off()


snv <- data.frame(cbind(snv_quartet[,c(1,3)],snv_ref[,c(3)]))
colnames(snv) <- c('caller','mendelian','f1')
snv$score <- snv$mendelian*snv$f1
snv_long <- melt(snv)
snv_long$sd <- c(snv_quartet$sd,snv_ref$sd)
p<-ggplot(snv_long, aes(x=caller, y=value, group=variable,color=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable))+
  xlim(c('sentieon','RTG','HC','Samtools','FreeBayes','SNVer','Varscan','ISAAC'))+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))
pdf('snv_mendelian_f1_lineplot.pdf',height=3,width=4)
plot(p)
dev.off()




snv <- data.frame(cbind(snv_quartet[,c(1,3)],snv_ref[,c(3)]))
colnames(snv) <- c('caller','mendelian','f1')
snv$score <- snv$mendelian*snv$f1
snv_long <- melt(snv)
snv_long$sd <- c(snv_quartet$sd,snv_ref$sd)
p<-ggplot(snv_long, aes(x=caller, y=value, group=variable,color=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable))+
  xlim(c('RTG','sentieon','HC','ISAAC','Varscan','FreeBayes','Samtools','SNVer'))+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(0.05))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))
pdf('indel_mendelian_f1_lineplot.pdf',height=3,width=4)
plot(p)
dev.off()



#sv
dat <- read.table('SV_precision_recall.txt',header=T)

ins <- dat[which(dat$SVTYPE == 'INS'),]
p <- ggplot(snv, aes(x=f1, y=mendelian,shape=mapper, color=caller)) + geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))

pdf('nctr_snv_scatter.pdf',height=4,width=5)
plot(p)
dev.off()

indel <- nctr[which(nctr$type == 'INDEL'),]
p <- ggplot(indel, aes(x=f1, y=mendelian,shape=mapper, color=caller)) + geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))
pdf('nctr_indel_scatter.pdf',height=4,width=5)
plot(p)
dev.off()



#
dat <- read.table('nctr_precision_recall.txt',header=T)
dat$f1 <- 2*dat$precision*dat$recall/(dat$precision+dat$recall)

ref <- data_summary(dat, varname="f1", 
                    groupnames=c("mapper",'caller','type'))
colnames(ref) <- c('mapper','caller','type','score','sd')
ref$tag <- 'F1'

dat2 <- read.table('nctr_quartet.txt',header=T)
quartet <- data_summary(dat2, varname="mendelian", 
                        groupnames=c("mapper",'caller','type'))

colnames(quartet) <- c('mapper','caller','type','score','sd')
quartet$tag <-'MCR'

df <- data.frame(rbind(ref,quartet))
df$score2 <- ifelse(df$tag == "MCR", -1 * df$score, df$score)

snv <- df[which(df$type=='INDEL'),]
snv$x <- paste(snv$caller,snv$mapper,sep="_")

p <- ggplot(data = snv) + geom_bar(aes(x=x,y=score2,fill=tag),stat="identity",position="identity") +
  scale_y_continuous(labels=abs)+
  geom_errorbar(aes(x=x,ymin=score2-sd, ymax=score2+sd), width=.2,
                position=position_dodge(0.05))+
  xlim(c('RTG_RTG','sentieon_sentieon','HC_Bowtie2','HC_BWA','HC_ISAAC','HC_stampy','ISAAC_Bowtie2','ISAAC_BWA','ISAAC_ISAAC','ISAAC_stampy',
         'Varscan_Bowtie2','Varscan_BWA','Varscan_ISAAC','Varscan_stampy','FreeBayes_Bowtie2','FreeBayes_BWA','FreeBayes_ISAAC','FreeBayes_stampy',
         'Samtools_Bowtie2','Samtools_BWA','Samtools_ISAAC','Samtools_stampy','SNVer_Bowtie2','SNVer_BWA','SNVer_ISAAC','SNVer_stampy'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))

pdf('indel_bidirectional_barplot.pdf',height=4,width=7)
plot(p)
dev.off()



###
# panel a the same pipeline but produced in different labs and technologies
pr_dat <- read.csv('Quartet_DNA_v202104/QC_metrics/snv.csv',header=T)
snv_ref <- data_summary(pr_dat, varname="METRIC.F1_Score", 
                    groupnames=c("platform",'library'))
colnames(snv_ref) <- c('platform','library','score','sd')
snv_ref$type <- 'SNV'

pr_dat <- read.csv('Quartet_DNA_v202104/QC_metrics/indel.csv',header=T)
indel_ref <- data_summary(pr_dat, varname="METRIC.F1_Score", 
                        groupnames=c("platform",'library'))
colnames(indel_ref) <- c('platform','library','score','sd')
indel_ref$type <- 'INDEL'
ref <- data.frame(rbind(snv_ref,indel_ref))
ref$tag <- 'F1'

dat <- read.table('quartet_all.txt',header = T)
quartet <- data_summary(dat, varname="mcr", 
                        groupnames=c("platform",'library','type'))

colnames(quartet) <- c('platform','library','type','score','sd')
quartet$tag <-'MCR'


ref <- ref[,c('platform','library','type','score','sd','tag')]
df <- data.frame(rbind(ref,quartet))
df$score2 <- ifelse(df$tag == "MCR", -1 * df$score, df$score)

snv <- df[which(df$type=='SNV'),]
snv$x <- paste(snv$library,snv$platform,sep="_")

p <- ggplot(data = snv) + geom_bar(aes(x=x,y=score2,fill=tag),stat="identity",position="identity") +
  ylim(c(-1,1))+
  geom_errorbar(aes(x=x,ymin=score2-sd, ymax=score2+sd), width=.2,
                position=position_dodge(0.05))+
#  xlim(c('RTG_RTG','sentieon_sentieon','HC_Bowtie2','HC_BWA','HC_ISAAC','HC_stampy','ISAAC_Bowtie2','ISAAC_BWA','ISAAC_ISAAC','ISAAC_stampy',
#         'Varscan_Bowtie2','Varscan_BWA','Varscan_ISAAC','Varscan_stampy','FreeBayes_Bowtie2','FreeBayes_BWA','FreeBayes_ISAAC','FreeBayes_stampy',
#         'Samtools_Bowtie2','Samtools_BWA','Samtools_ISAAC','Samtools_stampy','SNVer_Bowtie2','SNVer_BWA','SNVer_ISAAC','SNVer_stampy'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))


pdf('snv_platform_bidirectional_barplot.pdf',height=4,width=4)
plot(p)
dev.off()


####
#SV
dat <- read.delim('sv_quartet_metrics.txt',header=T)
dat <- dat[,-c(3,4,7)]
dat_long <- melt(dat)

#
dat_long$score2 <- ifelse(dat_long$variable == "mendelian", -1 * dat_long$value, dat_long$value)

dat_long$x <- paste(dat_long$Tech,dat_long$Pipelines,sep="_")

ins <- dat_long[which(dat_long$Type == 'INS'),]

p <- ggplot(data = ins) + geom_bar(aes(x=x,y=score2,fill=variable),stat="identity",position="identity") +
  ylim(c(-0.5,1))+
  #  xlim(c('RTG_RTG','sentieon_sentieon','HC_Bowtie2','HC_BWA','HC_ISAAC','HC_stampy','ISAAC_Bowtie2','ISAAC_BWA','ISAAC_ISAAC','ISAAC_stampy',
  #         'Varscan_Bowtie2','Varscan_BWA','Varscan_ISAAC','Varscan_stampy','FreeBayes_Bowtie2','FreeBayes_BWA','FreeBayes_ISAAC','FreeBayes_stampy',
  #         'Samtools_Bowtie2','Samtools_BWA','Samtools_ISAAC','Samtools_stampy','SNVer_Bowtie2','SNVer_BWA','SNVer_ISAAC','SNVer_stampy'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))


pdf('ins_long_reads_bidirectional_barplot.pdf',height=4.5,width=7)
plot(p)
dev.off()

######
all_value <- c(dat$F1_score *dat$mendelian)
dat_long$caller_order <- 1
dat_long$caller_order[which(dat_long$caller == 'nanosv')] <- 2
dat_long$caller_order[which(dat_long$caller == 'svim')] <- 3
dat_long$caller_order[which(dat_long$caller == 'sniffles')] <- 4
dat_long$caller_order[which(dat_long$caller == 'cutesv')] <- 5
dat_long$value_order <- c(all_value,all_value)
dat_long$plot_order <- paste(dat_long$Tech,dat_long$value_order,sep="_")
dat_long$variable <- as.character(dat_long$variable)

h1 <- ggplot(dat_long, aes(x=plot_order, y=value, group=variable, color=variable)) +
  geom_line()+
  geom_point() +
  scale_fill_brewer(palette="Paired")+
  facet_grid(~Type, scales = 'free_x', space = 'free_x')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size = 13,face = "bold",color = 'black'))

h4 <- ggplot(dat_long)+
  geom_bar(mapping = aes(x = plot_order, y = 1, fill = caller), 
           stat = "identity", 
           width = 1)+
  scale_fill_brewer(palette="Paired")+
  facet_grid(~Type, scales = 'free_x', space = 'free_x')+
  theme_void()+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing.x = unit(1, "mm"))


h3 <- ggplot(dat_long)+
  geom_bar(mapping = aes(x = plot_order, y = 1, fill = mapper), 
           stat = "identity", 
           width = 1)+
  facet_grid(~Type, scales = 'free_x', space = 'free_x')+
  theme_void()+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  scale_fill_brewer(palette="Set1")+
  theme(panel.spacing.x = unit(1, "mm"))

h2 <- ggplot(dat_long)+
  geom_bar(mapping = aes(x = plot_order, y = 1, fill = Tech), 
           stat = "identity", 
           width = 1)+
  facet_grid(~Type, scales = 'free_x', space = 'free_x')+
  theme_void()+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.spacing.x = unit(1, "mm"))


legend <- plot_grid(get_legend(h4), get_legend(h3), get_legend(h2), get_legend(h1), ncol = 1)
h1 <- h1 + theme(legend.position = "none")
h2 <- h2 + theme(legend.position = "none")
h3 <- h3 + theme(legend.position = "none")
h4 <- h4 + theme(legend.position = "none")
plot <- plot_grid(h1,h4,h3,h2,align = "v", ncol = 1, axis = "tb", rel_heights = c(15,0.4,0.4,0.4))
all_plot <- plot_grid(plot, legend, nrow = 1, rel_widths = c(10, 1.5))
plot(all_plot)
ggsave(all_plot,filename = "sv_metrics.pdf",width = 7,height = 3)

######## sv short-reads
dat <- read.table('sv_short_reads.txt',header=T)
dat <- dat[,c(2,3,4,10,16)]
dat_long <- melt(dat)

dat_long$score2 <- ifelse(dat_long$variable == "MCR", -1 * dat_long$value, dat_long$value)

dat_long$x <- paste(dat_long$Aligner,dat_long$Caller,sep="_")
ins <- dat_long[which(dat_long$SVTYPE == 'DEL'),]

p <- ggplot(data = ins) + geom_bar(aes(x=x,y=score2,fill=variable),stat="identity",position="identity") +
  ylim(c(-0.7,0.7))+
  #  xlim(c('RTG_RTG','sentieon_sentieon','HC_Bowtie2','HC_BWA','HC_ISAAC','HC_stampy','ISAAC_Bowtie2','ISAAC_BWA','ISAAC_ISAAC','ISAAC_stampy',
  #         'Varscan_Bowtie2','Varscan_BWA','Varscan_ISAAC','Varscan_stampy','FreeBayes_Bowtie2','FreeBayes_BWA','FreeBayes_ISAAC','FreeBayes_stampy',
  #         'Samtools_Bowtie2','Samtools_BWA','Samtools_ISAAC','Samtools_stampy','SNVer_Bowtie2','SNVer_BWA','SNVer_ISAAC','SNVer_stampy'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11,color='black'))+ # text on X axis
  theme(axis.text.y = element_text(size = 11,color='black'))+
  theme(axis.title.y = element_text(size = 14))

pdf('del_short_reads_bidirectional_barplot.pdf',height=4,width=5)
plot(p)
dev.off()

