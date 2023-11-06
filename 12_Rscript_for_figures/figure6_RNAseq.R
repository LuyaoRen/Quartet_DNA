library(ggplot2)
dat <- read.csv('targeted_all.csv')
p <- ggplot(dat, aes(x=METRIC.Precision, y=METRIC.Recall, color=caller, shape=Type)) +
  geom_point()+
  theme_light()+
  ylim(0.6,1)+
  xlim(0,1)+
  ylab('Recall')+
  xlab('Precision')

pdf('deepvariant_gatk_precision_recall.pdf',height = 2.5,width = 4)
plot(p)
dev.off()

dat <- dat[which(dat$Type == 'SNP'),c(4,7,8,19,20,21)]
dat$Error <- dat$QUERY.FP-dat$FP_rnaEditing
dat$name <- paste(dat$caller,dat$sample,dat$rep,sep="_")

dat <- dat[,c(1,3,7,8)]
colnames(dat) <- c('True Positive',"RNA editing","Potential artifacts",'sample')
library(reshape2)
dat_long <- melt(dat)
dat_long$variable <- factor(dat_long$variable, levels=c('Potential artifacts', 'RNA editing', 'True Positive'))
p <- ggplot(dat_long, aes(x=sample, y=value, fill=variable)) +
  geom_bar(position='stack', stat='identity')+
  theme_light()+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size = 8,color="black"))+ 
  theme(axis.text.y = element_text(size = 9))+
  theme(axis.title.y = element_text(size = 9))+
  ylab('SNV number')+
  xlab('')+
  scale_fill_manual(values=c('grey', '#22bb9c', 'pink'))
  
pdf('deepvariant_gatk_snv_number.pdf',height = 3,width = 5)
plot(p)
dev.off()