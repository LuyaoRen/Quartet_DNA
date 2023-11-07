
MNC=rbind(M5,N5,F5)
MNC$depth=MNC$V6+MNC$V7
MNC$AF=MNC$V7/MNC$depth


SVF=read.delim("clipboard",header=F,sep = "\t")
SVF$class="Filter"
SVN=read.delim("clipboard",header=F,sep = "\t")
SVN$class="NM"
SVM=read.delim("clipboard",header=F,sep = "\t")
SVM$class="M"
FMN=rbind(SVF,SVN,SVM)
FMN$AF=FMN$V7/(FMN$V6+FMN$V7)
FMN$depth=FMN$V6+FMN$V7
FMN$class=factor(FMN$class,levels=c("Filter","NM","M"),labels=c("Non-repeatable","Technology or pipeline \nconcordance","Mendelian concordance"))
col=c("#DF8F44FF","#868686FF","#00A1D5FF")
#AF
ggplot(FMN, aes(x=AF, fill=class)) +
  geom_density(alpha=0.5)+scale_fill_manual(name="SV class",values=col) + theme_classic() + #theme_bw() +
  xlab("AF") + ylab("Density") +
  theme(text=element_text(color='black'),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))
#PDF 5/3.5
#RE
ggplot(FMN[FMN$depth<=100,], aes(x=V7, fill=class)) +
  geom_density(alpha=0.5)+scale_fill_manual(name="SV class",values=col) + theme_classic() + #theme_bw() +
  xlab("Reads supporting SV") + ylab("Density") +
  theme(text=element_text(color='black'),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12))
#PDF 5/3.5

#Depth
ggplot(FMN[FMN$depth<=100,], aes(x=depth, fill=class)) +
  geom_density(alpha=0.5)+scale_fill_manual(values=col)
