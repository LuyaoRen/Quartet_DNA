library(ggplot2)
library(reshape2)

###boxplot
snv <- read.delim('performance_metrics.txt',header=T)
sv <- read.delim('performance_metrics_SV_GT_F1.txt',header=T)

snv <- snv[,c(3,4,6,16,15)]
sv <- sv[,c(3,4,6,10,11)]

colnames(snv) <- c('library','platform','Type','Precision','Recall')

dat <- data.frame(rbind(snv,sv))
dat$cate <- paste(dat$library,dat$platform,sep="_")
dat$Type_f = factor(dat$Type, levels=c('SNP','INDEL','DEL','INS'))
p <- ggplot(dat, aes(x = cate, y = Recall, color=cate))+
  geom_boxplot()+
  facet_wrap(~ Type_f, nrow=1)+
  ylab('Recall')+
  xlab('')+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  theme(axis.text.x = element_blank())+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size=8, color="black"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_y_log10()
  
pdf('recall.pdf',height=2,width=10)  
plot(p)
dev.off()

#mendelian

snv <- read.delim('performance_metrics.txt',header=T)
sv <- read.delim('SV_noGT_MCR.txt',header=T)

snv <- snv[,c(3,4,6,43,44)]
sv <- sv[,c(3,4,5,22,21)]

colnames(snv) <- c('library','platform','Type','Quartet_inside','Quartet_outside')
snv$Quartet_inside <- 1-snv$Quartet_inside
snv$Quartet_outside <- 1-snv$Quartet_outside
colnames(sv) <- c('library','platform','Type','Quartet_inside','Quartet_outside')


dat <- data.frame(rbind(snv,sv))
dat$cate <- paste(dat$library,dat$platform,sep="_")
dat$Type_f = factor(dat$Type, levels=c('SNP','INDEL','DEL','INS'))
p <- ggplot(dat, aes(x = cate, y = Quartet_inside, color=cate))+
  geom_boxplot()+
  facet_wrap(~ Type_f, nrow=1,scales='free_x')+
  ylab('Menlian violation rate')+
  xlab('')+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  theme(axis.text.x = element_blank())+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size=8, color="black"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

pdf('quartet_inside.pdf',height=2,width=10)  
plot(p)
dev.off()

dat <- dat[,c(7,4,5)]
dat_long <- melt(dat)
dat_long$cate <- paste(dat_long$Type_f,dat_long$variable,sep="_")
p <- ggplot(dat_long, aes(x = cate, y = value, color=variable))+
  geom_boxplot()+
  #facet_wrap(~ Type_f, nrow=1,scales='free_x')+
  ylab('Menlian violation rate')+
  xlab('')+
  xlim(c('SNP_Quartet_inside','SNP_Quartet_outside','INDEL_Quartet_inside','INDEL_Quartet_outside','DEL_Quartet_inside','DEL_Quartet_outside','INS_Quartet_inside','INS_Quartet_outside'))+
  theme_bw()+
  scale_color_brewer(palette="Paired")+
  theme(axis.text.x = element_blank())+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size=8, color="black"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

pdf('inside_outside_boxplot.pdf',height=2.5,width=5)  
plot(p)
dev.off()


#library
dat <- dat[,c(1,7,6,4,5)]
dat_long <- melt(dat)
dat_long$cate <- paste(dat_long$Type_f,dat_long$variable,sep="_")

p <- ggplot(dat_long, aes(x = variable, y = value, color=variable))+
  geom_boxplot()+
  facet_wrap(cate ~ Type_f, nrow = 7,scales='free_x')+
  ylab('Menlian violation rate')+
  xlab('')+
  #xlim(c('SNP_Quartet_inside','SNP_Quartet_outside','INDEL_Quartet_inside','INDEL_Quartet_outside','DEL_Quartet_inside','DEL_Quartet_outside','INS_Quartet_inside','INS_Quartet_outside'))+
  theme_bw()+
  scale_color_manual(values = c('#f6a945','#563405'))+
  #theme(axis.text.x = element_blank())+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))+
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size=8, color="black"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
pdf('library_platform_region.pdf',height = 10,width=8)
plot(p)
dev.off()

##inside vs outside

snv <- read.delim('performance_metrics.txt',header=T)
sv <- read.delim('SV_noGT_MCR.txt',header=T)

snv <- snv[,c(3,4,6,16,43)]
sv <- sv[,c(3,4,5,25,22)]

colnames(snv) <- c('library','platform','Type','precision','mcr')
#snv$Quartet_inside <- 1-snv$Quartet_inside
#snv$Quartet_outside <- 1-snv$Quartet_outside
colnames(sv) <- c('library','platform','Type','precision','mcr')
sv$mcr <- 1-sv$mcr
dat <- data.frame(rbind(snv,sv))
dat$cate <- paste(dat$library,dat$platform,sep="_")
dat$Type_f = factor(dat$Type, levels=c('SNP','INDEL','DEL','INS'))


p <- ggplot(dat, aes(x=precision, y=mcr,color=platform,shape=library)) + geom_point()+
  facet_wrap(~ Type_f, nrow=1,scales='free')+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  theme(axis.text.x = element_text(size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))+
  theme(strip.text.x = element_text(size=8, color="black"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf('inside_outside_scatter_plot.pdf',height=2,width=10)  
plot(p)
dev.off()













#### pedigree
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

sub <- dat[,c(1,2,3,4,5,6,34,35,36,37,38,39,40,41)]
#twins
sub_dat_summary <- data_summary(sub, varname="twins_outside_consistency.1", 
             groupnames=c("batch",'library','platform','Type'))
indel <- sub_dat_summary[which(sub_dat_summary$Type=='SNP'),]


p <- ggplot(indel, aes(x=batch, y=twins_outside_consistency.1*1000000,fill=platform)) + 
  geom_bar(stat="identity", color='black', position=position_dodge()) +
  geom_errorbar(aes(ymin=twins_outside_consistency.1*1000000-sd*1000000, ymax=twins_outside_consistency.1*1000000+sd*1000000), width=.2,
                position=position_dodge(.9))+
  scale_fill_brewer(palette="Dark2")+
  ylab('False Positive per Mb')+
  xlab('')+
  ylim(c(0,480))+
  xlim(c('BGI_SEQ500_BGI_1','BGI_SEQ500_BGI_2','BGI_SEQ500_BGI_3','ILM_Nova_BRG','ILM_Nova_GAC','ILM_Nova_NVG','ILM_Nova_WUX','ILM_XTen_WUX_4','ILM_XTen_WUX_5','ILM_XTen_WUX_6',
         'BGI_T5_SZ_1','BGI_T5_SZ_2','BGI_T5_SZ_3','MGI_T7RS_MingMa','MGI_T7RS_WeGene','illumina_Nova_MingMa','illumina_Novaseq6000_KingMed'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))

pdf('snv_twin_outside.pdf',height=3.5,width=5.5)
plot(p)
dev.off()


#trio outside
sub <- dat[,c(2,3,4,5,6,39,40)]
sub_long <- melt(sub)
sub_long <- na.omit(sub_long)
sub_dat_summary <- data_summary(sub_long, varname="value", 
                                groupnames=c("batch",'library','platform','Type'))

indel <- sub_dat_summary[which(sub_dat_summary$Type=='INDEL'),]


p <- ggplot(indel, aes(x=batch, y=value*1000000,fill=platform)) + 
  geom_bar(stat="identity", color='black', position=position_dodge()) +
  geom_errorbar(aes(ymin=value*1000000-sd*1000000, ymax=value*1000000+sd*1000000), width=.2,
                position=position_dodge(.9))+
  scale_fill_brewer(palette="Dark2")+
  ylab('False Positive per Mb')+
  xlab('')+
  ylim(c(0,480))+
  xlim(c('BGI_SEQ500_BGI_1','BGI_SEQ500_BGI_2','BGI_SEQ500_BGI_3','ILM_Nova_BRG','ILM_Nova_GAC','ILM_Nova_NVG','ILM_Nova_WUX','ILM_XTen_WUX_4','ILM_XTen_WUX_5','ILM_XTen_WUX_6',
         'BGI_T5_SZ_1','BGI_T5_SZ_2','BGI_T5_SZ_3','MGI_T7RS_MingMa','MGI_T7RS_WeGene','illumina_Nova_MingMa','illumina_Novaseq6000_KingMed'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))

pdf('indel_trio_outside.pdf',height=3.5,width=5.5)
plot(p)
dev.off()

#quartet outside


sub <- dat[,c(1,2,3,4,5,6,41)]
#twins
sub_dat_summary <- data_summary(sub, varname="quartet_outside.1", 
                                groupnames=c("batch",'library','platform','Type'))
indel <- sub_dat_summary[which(sub_dat_summary$Type=='SNP'),]


p <- ggplot(indel, aes(x=batch, y=quartet_outside.1*1000000,fill=platform)) + 
  geom_bar(stat="identity", color='black', position=position_dodge()) +
  geom_errorbar(aes(ymin=quartet_outside.1*1000000-sd*1000000, ymax=quartet_outside.1*1000000+sd*1000000), width=.2,
                position=position_dodge(.9))+
  scale_fill_brewer(palette="Dark2")+
  ylab('False Positive per Mb')+
  xlab('')+
  xlim(c('BGI_SEQ500_BGI_1','BGI_SEQ500_BGI_2','BGI_SEQ500_BGI_3','ILM_Nova_BRG','ILM_Nova_GAC','ILM_Nova_NVG','ILM_Nova_WUX','ILM_XTen_WUX_4','ILM_XTen_WUX_5','ILM_XTen_WUX_6',
         'BGI_T5_SZ_1','BGI_T5_SZ_2','BGI_T5_SZ_3','MGI_T7RS_MingMa','MGI_T7RS_WeGene','illumina_Nova_MingMa','illumina_Novaseq6000_KingMed'))+
  ylim(c(0,480))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))

pdf('snv_quartet_outside.pdf',height=3.5,width=5.5)
plot(p)
dev.off()


###sv
dat <- read.table('SV_noGT_MCR.txt',header=T)
#twins
ins <- dat[which(dat$Type=='INS'),]


p <- ggplot(ins, aes(x=batch, y=D5_D6_FPR,fill=platform)) + 
  geom_bar(stat="identity", color='black', position=position_dodge()) +
  scale_fill_brewer(palette="Dark2")+
  ylab('False Positive per Mb')+
  xlab('')+
  ylim(c(0,1.3))+
  xlim(c('BGI_SEQ500_BGI_1','BGI_SEQ500_BGI_2','BGI_SEQ500_BGI_3','ILM_Nova_BRG','ILM_Nova_GAC','ILM_Nova_NVG','ILM_Nova_WUX','ILM_XTen_WUX_4','ILM_XTen_WUX_5','ILM_XTen_WUX_6',
         'BGI_T5_SZ_1','BGI_T5_SZ_2','BGI_T5_SZ_3','MGI_T7RS_MingMa','MGI_T7RS_WeGene','illumina_Nova_MingMa','illumina_Novaseq6000_KingMed'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))

pdf('sv_ins_twin_outside.pdf',height=3.5,width=5.5)
plot(p)
dev.off()


#trio outside

p <- ggplot(ins, aes(x=batch, y=D5_F7_M8_FPR,fill=platform)) + 
  geom_bar(stat="identity", color='black', position=position_dodge()) +
  scale_fill_brewer(palette="Dark2")+
  ylab('False Positive per Mb')+
  xlab('')+
  ylim(c(0,1.3))+
  xlim(c('BGI_SEQ500_BGI_1','BGI_SEQ500_BGI_2','BGI_SEQ500_BGI_3','ILM_Nova_BRG','ILM_Nova_GAC','ILM_Nova_NVG','ILM_Nova_WUX','ILM_XTen_WUX_4','ILM_XTen_WUX_5','ILM_XTen_WUX_6',
         'BGI_T5_SZ_1','BGI_T5_SZ_2','BGI_T5_SZ_3','MGI_T7RS_MingMa','MGI_T7RS_WeGene','illumina_Nova_MingMa','illumina_Novaseq6000_KingMed'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))

pdf('sv_ins_trio_outside.pdf',height=3.5,width=5.5)
plot(p)
dev.off()

#quartet outside

p <- ggplot(ins, aes(x=batch, y=Quartet_FPR,fill=platform)) + 
  geom_bar(stat="identity", color='black', position=position_dodge()) +
  scale_fill_brewer(palette="Dark2")+
  ylab('False Positive per Mb')+
  xlab('')+
  xlim(c('BGI_SEQ500_BGI_1','BGI_SEQ500_BGI_2','BGI_SEQ500_BGI_3','ILM_Nova_BRG','ILM_Nova_GAC','ILM_Nova_NVG','ILM_Nova_WUX','ILM_XTen_WUX_4','ILM_XTen_WUX_5','ILM_XTen_WUX_6',
         'BGI_T5_SZ_1','BGI_T5_SZ_2','BGI_T5_SZ_3','MGI_T7RS_MingMa','MGI_T7RS_WeGene','illumina_Nova_MingMa','illumina_Novaseq6000_KingMed'))+
  ylim(c(0,1.3))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 11))+ # text on X axis
  theme(axis.text.y = element_text(size = 11))+
  theme(axis.title.y = element_text(size = 14))

pdf('sv_ins_quartet_outside.pdf',height=3.5,width=5.5)
plot(p)
dev.off()





