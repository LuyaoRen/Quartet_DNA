library(ggplot2)
library(RColorBrewer)
display.brewer.all()

a=read.delim("clipboard",header=F)
b=a[,c(5,6)]
c=b
c$V6=abs(b$V6)
LCL5_LCL6_f_nosame=c
colnames(LCL5_LCL6_f_nosame)=c("V1","V2")
LCL5_LCL6_f_nosame_total=LCL5_LCL6_f_nosame #5SV
col=brewer.pal(8,"Set1")[c(1,3,4)]
colnames(LCL5_LCL6_f_nosame_total)=c("V1","V2")
LCL5_LCL6_f_nosame_1K=LCL5_LCL6_f_nosame_total[LCL5_LCL6_f_nosame_total$V2<1000 & LCL5_LCL6_f_nosame_total$V2>=50,]
ggplot(data = LCL5_LCL6_f_nosame_1K,aes(x=V2,fill=V1))+geom_histogram(bins=50)+
  ggtitle("")+labs(x="SV length(<1kb)")+theme_bw()+theme(text=element_text(color='black', size=18),legend.position=c(0.85,0.8))+
  guides(fill = guide_legend(title = NULL))+
  theme(panel.grid=element_line(colour="white"))+
  scale_fill_manual(values=c("#0066CC", "#FF9966"))

LCL5_LCL6_f_nosame_1K_10K=LCL5_LCL6_f_nosame_total[10000>LCL5_LCL6_f_nosame_total$V2 & LCL5_LCL6_f_nosame_total$V2>=1000,]
ggplot(data=LCL5_LCL6_f_nosame_1K_10K,aes(x=V2,fill=V1))+geom_histogram(bins=50)+
  ggtitle("")+labs(x="SV length(1kb~10kb)")+theme_bw()+theme(text=element_text(color='black', size=18),legend.position=c(0.85,0.8))+
  guides(fill = guide_legend(title = NULL))+ylim(0,200)+
  theme(panel.grid=element_line(colour="white"))+
  scale_fill_manual(values=c("#0066CC", "#FF9966"))

LCL5_LCL6_f_nosame_1K$tag=rep("1K",nrow(LCL5_LCL6_f_nosame_1K))
LCL5_LCL6_f_nosame_1K_10K$tag=rep("1K_10K",nrow(LCL5_LCL6_f_nosame_1K_10K))
M=rbind(LCL5_LCL6_f_nosame_1K,LCL5_LCL6_f_nosame_1K_10K)
ggplot(data=M,aes(x=V2,fill=V1))+geom_histogram(bins=50)+facet_wrap(tag~.,scales = "free")+
  ggtitle("")+theme_bw()+theme(text=element_text(color='black', size=18),legend.position=c(0.85,0.8))+
  guides(fill = guide_legend(title = NULL))+
  theme(panel.grid=element_line(colour="white"))+
  scale_fill_manual(values=c("#0066CC", "#FF9966"))



