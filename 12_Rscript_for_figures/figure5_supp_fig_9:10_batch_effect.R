setwd('Quartet_DNA_v202112/batch_effect/')
library(ggplot2)
dat_inside <- read.delim('merged.inside.chr1.indel.number.txt',header=F)
dat_outside <- read.delim('merged.outside.chr1.indel.number.txt',header=F)
#dat <- data.frame(rbind(dat_inside,dat_outside))
#dat2 <- dat[,-c(1,2,3,4)]

dat2 <- dat_inside[,-c(1,2,3,4)]
res.pca <- prcomp(t(dat2))
df_out <- as.data.frame(res.pca$x)
df_out$batch <- c('BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI',
                  'BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ',
                  'illumina_Nova_MingMa','illumina_Nova_MingMa','illumina_Nova_MingMa','illumina_Nova_MingMa',
                  'illumina_Novaseq6000_KingMed','illumina_Novaseq6000_KingMed','illumina_Novaseq6000_KingMed','illumina_Novaseq6000_KingMed',
                  'ILM_Nova_BRG','ILM_Nova_BRG','ILM_Nova_BRG','ILM_Nova_BRG',
                  'ILM_Nova_GAC','ILM_Nova_GAC','ILM_Nova_GAC','ILM_Nova_GAC',
                  'ILM_Nova_NVG','ILM_Nova_NVG','ILM_Nova_NVG','ILM_Nova_NVG',
                  'ILM_Nova_WUX','ILM_Nova_WUX','ILM_Nova_WUX','ILM_Nova_WUX',
                  'ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX',
                  'MGI_T7RS_MingMa','MGI_T7RS_MingMa','MGI_T7RS_MingMa','MGI_T7RS_MingMa',
                  'MGI_T7RS_WeGene','MGI_T7RS_WeGene','MGI_T7RS_WeGene','MGI_T7RS_WeGene'
                  )
df_out$sample <- c('LCL5','LCL5','LCL5','LCL6','LCL6','LCL6','LCL7','LCL7','LCL7','LCL8','LCL8','LCL8',
                   'LCL5','LCL5','LCL5','LCL6','LCL6','LCL6','LCL7','LCL7','LCL7','LCL8','LCL8','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL5','LCL5','LCL6','LCL6','LCL6','LCL7','LCL7','LCL7','LCL8','LCL8','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8')
pca_Variance<-summary(res.pca)$importance[2,]
p <- ggplot(df_out,aes(x=PC2,y=PC3,color=sample,shape=batch))+
  geom_point()+
  scale_color_manual(values = c('#4CC3D9','#7BC8A4','#FFC65D','#F16745'))+
  scale_shape_manual(values=c(15, 16, 17,18,0,1,2,3,4,7,8,9))+
  xlab(paste("PC2  ", round(pca_Variance[2] * 100,2),"%",sep=""))+
  ylab(paste("PC3  ", round(pca_Variance[3] * 100,2),"%",sep=""))+
  theme_bw()
pdf('indel.pc23.inside.pdf',height=3,width=5.5)
plot(p)
dev.off()

##sv
##insertion 
dat <- read.delim('sv/sv.ins.all.txt',header=F)
dat2 <- dat[,-c(1,2,3)]
res.pca <- prcomp(t(dat2))
df_out <- as.data.frame(res.pca$x)
df_out$batch <- c('BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI','BGI_SEQ500_BGI',
                  'BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ','BGI_T5_SZ',
                  'illumina_Nova_MingMa','illumina_Nova_MingMa','illumina_Nova_MingMa','illumina_Nova_MingMa',
                  'illumina_Novaseq6000_KingMed','illumina_Novaseq6000_KingMed','illumina_Novaseq6000_KingMed','illumina_Novaseq6000_KingMed',
                  'ILM_Nova_BRG','ILM_Nova_BRG','ILM_Nova_BRG','ILM_Nova_BRG',
                  'ILM_Nova_GAC','ILM_Nova_GAC','ILM_Nova_GAC','ILM_Nova_GAC',
                  'ILM_Nova_NVG','ILM_Nova_NVG','ILM_Nova_NVG','ILM_Nova_NVG',
                  'ILM_Nova_WUX','ILM_Nova_WUX','ILM_Nova_WUX','ILM_Nova_WUX',
                  'ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX','ILM_XTen_WUX',
                  'MGI_T7RS_MingMa','MGI_T7RS_MingMa','MGI_T7RS_MingMa','MGI_T7RS_MingMa',
                  'MGI_T7RS_WeGene','MGI_T7RS_WeGene','MGI_T7RS_WeGene','MGI_T7RS_WeGene'
)
df_out$sample <- c('LCL5','LCL5','LCL5','LCL6','LCL6','LCL6','LCL7','LCL7','LCL7','LCL8','LCL8','LCL8',
                   'LCL5','LCL5','LCL5','LCL6','LCL6','LCL6','LCL7','LCL7','LCL7','LCL8','LCL8','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL5','LCL5','LCL6','LCL6','LCL6','LCL7','LCL7','LCL7','LCL8','LCL8','LCL8',
                   'LCL5','LCL6','LCL7','LCL8',
                   'LCL5','LCL6','LCL7','LCL8')
pca_Variance<-summary(res.pca)$importance[2,]
p <- ggplot(df_out,aes(x=PC1,y=PC2,color=sample,shape=batch))+
  geom_point()+
  scale_color_manual(values = c('#4CC3D9','#7BC8A4','#FFC65D','#F16745'))+
  scale_shape_manual(values=c(15, 16, 17,18,0,1,2,3,4,7,8,9))+
  xlab(paste("PC1  ", round(pca_Variance[1] * 100,2),"%",sep=""))+
  ylab(paste("PC2  ", round(pca_Variance[2] * 100,2),"%",sep=""))+
  theme_bw()
pdf('ins.pc12.pdf',height=3,width=5.5)
plot(p)
dev.off()

plot_ly(
  df_out, x = ~PC2, y = ~PC3, z = ~PC4, 
  color = ~batch
) %>%
  add_markers() %>%
  layout(
    scene = list(xaxis = list(title = 'PC2'),
                 yaxis = list(title = 'PC3'),
                 zaxis = list(title = 'PC4'))
  )

##### reproducibility
##snv
file_1 <- c('snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_2_20180328_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_3_20180328_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_2_20180328_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_3_20180328_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_2_20200716_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_3_20200716_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_2_20180328_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_3_20180328_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_2_20180328_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_3_20180328_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_2_20200716_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_3_20200716_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_2_20180328_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_3_20180328_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_2_20180328_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_3_20180328_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_2_20200716_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_3_20200716_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_2_20180328_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_3_20180328_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_2_20180328_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_3_20180328_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_2_20200716_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_3_20200716_hc.vcf.info.txt.all.snv')

file_2 <- c('snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_2_20200716_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_3_20200716_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_5_20180703_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_5_20180703_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_2_20200716_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_3_20200716_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_5_20180703_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_5_20180703_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_2_20200716_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_3_20200716_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_5_20180703_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_5_20180703_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_2_20200716_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_3_20200716_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_5_20180703_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_5_20180703_hc.vcf.info.txt.all.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_6_20180703_hc.vcf.info.txt.all.snv'
)

filtered_snv <- c()
for(i in c(1:24)){
  filtered_1_snv <- read.table(file_1[i],header=F)
  filtered_1_snv <- paste(filtered_1_snv$V1,filtered_1_snv$V2,filtered_1_snv$V3,filtered_1_snv$V4,filtered_1_snv$V5,sep="_")
  filtered_2_snv <- read.table(file_2[i],header=F)
  filtered_2_snv <- paste(filtered_2_snv$V1,filtered_2_snv$V2,filtered_2_snv$V3,filtered_2_snv$V4,filtered_2_snv$V5,sep="_")
  inter_filtered_snv <- intersect(filtered_1_snv,filtered_2_snv)
  union_filtered_snv <- union(filtered_1_snv,filtered_2_snv)
  R <- 1/2*((length(inter_filtered_snv)/length(filtered_1_snv)) + (length(inter_filtered_snv)/length(filtered_2_snv)))
  filtered_snv <- c(filtered_snv,R)
}


file_1 <- c('snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_2_20180328_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_3_20180328_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_2_20180328_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_3_20180328_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_2_20200716_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_3_20200716_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_2_20180328_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_3_20180328_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_2_20180328_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_3_20180328_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_2_20200716_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_3_20200716_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_2_20180328_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_3_20180328_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_2_20180328_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_3_20180328_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_2_20200716_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_3_20200716_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_2_20180328_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_3_20180328_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_2_20180328_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_3_20180328_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_2_20200716_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_3_20200716_hc.vcf.info.txt.all.indel')

file_2 <- c('snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_2_20200716_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_3_20200716_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_5_20180703_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_5_20180703_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_2_20200716_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_3_20200716_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_5_20180703_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_5_20180703_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_2_20200716_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_3_20200716_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_5_20180703_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_5_20180703_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_2_20200716_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_3_20200716_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_5_20180703_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_5_20180703_hc.vcf.info.txt.all.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_6_20180703_hc.vcf.info.txt.all.indel'
)
filtered_indel <- c()
for(i in c(1:24)){
  filtered_1_indel <- read.table(file_1[i],header=F)
  filtered_1_indel <- paste(filtered_1_indel$V1,filtered_1_indel$V2,filtered_1_indel$V3,filtered_1_indel$V4,filtered_1_indel$V5,sep="_")
  filtered_2_indel <- read.table(file_2[i],header=F)
  filtered_2_indel <- paste(filtered_2_indel$V1,filtered_2_indel$V2,filtered_2_indel$V3,filtered_2_indel$V4,filtered_2_indel$V5,sep="_")
  inter_filtered_indel <- intersect(filtered_1_indel,filtered_2_indel)
  union_filtered_indel <- union(filtered_1_indel,filtered_2_indel)
  R <- 1/2*((length(inter_filtered_indel)/length(filtered_1_indel)) + (length(inter_filtered_indel)/length(filtered_2_indel)))
  filtered_indel <- c(filtered_indel,R)
}


###all

file_1 <- c('snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_2_20180328_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_3_20180328_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_2_20180328_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_3_20180328_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_2_20200716_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_3_20200716_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_2_20180328_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_3_20180328_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_2_20180328_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_3_20180328_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_2_20200716_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_3_20200716_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_2_20180328_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_3_20180328_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_2_20180328_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_3_20180328_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_2_20200716_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_3_20200716_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_2_20180328_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_3_20180328_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_2_20180328_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_3_20180328_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_2_20200716_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_3_20200716_hc.vcf.info.txt.filtered.snv'
)

file_2 <- c('snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_2_20200716_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_3_20200716_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_5_20180703_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_5_20180703_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_2_20200716_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_3_20200716_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_5_20180703_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_5_20180703_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_2_20200716_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_3_20200716_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_5_20180703_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_5_20180703_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_2_20200716_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_3_20200716_hc.vcf.info.txt.filtered.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_5_20180703_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_6_20180703_hc.vcf.info.txt.all.snv',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_5_20180703_hc.vcf.info.txt.filtered.snv','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_6_20180703_hc.vcf.info.txt.all.snv'
            )
all_snv <- c()
for(i in c(1:24)){
  filtered_1_snv <- read.table(file_1[i],header=F)
  filtered_1_snv <- paste(filtered_1_snv$V1,filtered_1_snv$V2,filtered_1_snv$V3,filtered_1_snv$V4,filtered_1_snv$V5,sep="_")
  filtered_2_snv <- read.table(file_2[i],header=F)
  filtered_2_snv <- paste(filtered_2_snv$V1,filtered_2_snv$V2,filtered_2_snv$V3,filtered_2_snv$V4,filtered_2_snv$V5,sep="_")
  inter_filtered_snv <- intersect(filtered_1_snv,filtered_2_snv)
  union_filtered_snv <- union(filtered_1_snv,filtered_2_snv)
  R <- 1/2*((length(inter_filtered_snv)/length(filtered_1_snv)) + (length(inter_filtered_snv)/length(filtered_2_snv)))
  all_snv <- c(all_snv,R)
}


##all indel
file_1 <- c('snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_2_20180328_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_3_20180328_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_2_20180328_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL5_3_20180328_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_2_20200716_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_3_20200716_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_2_20180328_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_3_20180328_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_2_20180328_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL6_3_20180328_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_2_20200716_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_3_20200716_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_2_20180328_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_3_20180328_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_2_20180328_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL7_3_20180328_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_2_20200716_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_3_20200716_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_2_20180328_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_3_20180328_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_2_20180328_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_SEQ500_BGI_LCL8_3_20180328_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_2_20200716_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_3_20200716_hc.vcf.info.txt.filtered.indel'
)

file_2 <- c('snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_2_20200716_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL5_3_20200716_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_5_20180703_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_5_20180703_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL5_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_2_20200716_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL6_3_20200716_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_5_20180703_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_5_20180703_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL6_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_2_20200716_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL7_3_20200716_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_5_20180703_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_5_20180703_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL7_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_2_20200716_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_BGI_T5_SZ_LCL8_3_20200716_hc.vcf.info.txt.filtered.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_5_20180703_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_6_20180703_hc.vcf.info.txt.all.indel',
            'snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_5_20180703_hc.vcf.info.txt.filtered.indel','snv_pos/Quartet_DNA_ILM_XTen_WUX_LCL8_6_20180703_hc.vcf.info.txt.all.indel'
)
all_indel <- c()
for(i in c(1:24)){
  filtered_1_snv <- read.table(file_1[i],header=F)
  filtered_1_snv <- paste(filtered_1_snv$V1,filtered_1_snv$V2,filtered_1_snv$V3,filtered_1_snv$V4,filtered_1_snv$V5,sep="_")
  filtered_2_snv <- read.table(file_2[i],header=F)
  filtered_2_snv <- paste(filtered_2_snv$V1,filtered_2_snv$V2,filtered_2_snv$V3,filtered_2_snv$V4,filtered_2_snv$V5,sep="_")
  inter_filtered_snv <- intersect(filtered_1_snv,filtered_2_snv)
  union_filtered_snv <- union(filtered_1_snv,filtered_2_snv)
  R <- 1/2*((length(inter_filtered_snv)/length(filtered_1_snv)) + (length(inter_filtered_snv)/length(filtered_2_snv)))
  all_indel <- c(all_indel,R)
}

filtered_indel <- data.frame(filtered_indel)
filtered_indel$tag <- 'raw'
filtered_indel$type <- 'INDEL'
colnames(filtered_indel)[1] <- 'number'

filtered_snv <- data.frame(filtered_snv)
filtered_snv$tag <- 'raw'
filtered_snv$type <- 'SNV'
colnames(filtered_snv)[1] <- 'number'


all_indel <- data.frame(all_indel)
all_indel$tag <- 'filtered'
all_indel$type <- 'INDEL'
colnames(all_indel)[1] <- 'number'


all_snv <- data.frame(all_snv)
all_snv$tag <- 'filtered'
all_snv$type <- 'SNV'
colnames(all_snv)[1] <- 'number'


dat <- data.frame(rbind(filtered_indel,filtered_snv,all_indel,all_snv))
write.table(dat,'snv_reproducibility.filter.txt',sep="\t",col.names=F,row.names=F)
dat <- read.table('reproducibility.filter.txt',header=F)
dat$V3_f <- factor(dat$V3,levels=c('SNV','INDEL','DEL','INS'))
library(ggplot2)
p <- ggplot(dat, aes(x=V2, y=V1, fill=V2)) +
  geom_boxplot()+
  facet_wrap(~ V3_f, nrow=1)+
  xlim(c('raw','filtered'))+
  ylab('Reproducibility')+
  xlab('')+
  theme_light()+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11))
pdf('reproducibility.pdf',height = 2.5,width = 4)
plot(p)
dev.off()

###precision
dat <- read.delim('vcf_filter.txt',header=T)

dat <- dat[,c(4,5,8)]
dat_long <- melt(dat)
dat_long$Type_f <- factor(dat_long$Type,levels=c('SNP','INDEL','DEL','INS'))
p <- ggplot(dat_long, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot()+
  facet_wrap(~ Type_f, nrow=1)+
  xlim(c('before_precision','after_precision'))+
  ylab('Precision')+
  xlab('')+
  theme_light()+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11))
pdf('precison.pdf',height = 2.5,width = 5)
plot(p)
dev.off()


dat <- read.delim('vcf_filter.txt',header=T)
dat <- dat[,c(4,6,9)]
dat_long <- melt(dat)
dat_long$Type_f <- factor(dat_long$Type,levels=c('SNP','INDEL','DEL','INS'))
p <- ggplot(dat_long, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot()+
  facet_wrap(~ Type_f, nrow=1)+
  xlim(c('before_recall','after_recall'))+
  ylab('Recall')+
  xlab('')+
  theme_light()+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11))
pdf('recall.pdf',height = 2.5,width = 5)
plot(p)
dev.off()

####SV threshold

line_dat <-function(all_table,mc_table,tag,mc_number,mv_number,type){
  bgi_500_all <- read.table(all_table,header=F)
  bgi_500_mc <- read.table(mc_table,header=F)
  bgi_500_mc$tag <- 'MC'
  bgi_500 <- merge(bgi_500_all,bgi_500_mc,by=c('V1','V2'),all = TRUE)
  bgi_500 <- bgi_500[which(is.na(bgi_500$V3)==FALSE),]
  bgi_500$tag[which(is.na(bgi_500$tag)==TRUE)] <- 'MV'
  bgi_500 <- bgi_500[which(bgi_500$V3 == type),]
  
  dat_point <- c(seq(1,1000,5),997,998)
  result_qd <- matrix(NA,ncol=2,nrow = 202)
  
  for(i in 1:length(dat_point)){
    print(i)
    loc <- which(bgi_500$V4 > dat_point[i])
    sub_men <- bgi_500[loc,]
    result_qd[i,1] <- length(which(sub_men[,'tag'] == 'MC'))
    result_qd[i,2] <- length(which(sub_men[,'tag'] == 'MV'))
  }
  
  result_qd[,1] <- result_qd[,1]/mc_number
  result_qd[,2] <- result_qd[,2]/mv_number
  colnames(result_qd) <- c('mendelian_consistent','mendelian_violation')
  result_qd <- data.frame(result_qd)
  result_qd$site <- tag
  result_qd$para <- dat_point
  return(result_qd)
}


bgi_500 <- line_dat('SV/Quartet_DNA_BGI_SEQ500_BGI_LCL5_1_20180328_manta_all.txt','SV/Quartet_DNA_BGI_SEQ500_BGI_LCL5_1_M.pos.txt','BGI_SEQ500_BGI',3189,607,'MantaDEL')
bgi_t5 <- line_dat('SV/Quartet_DNA_BGI_T5_SZ_LCL5_1_20200716_manta_all.txt','SV/Quartet_DNA_BGI_T5_SZ_LCL5_1_M.pos.txt','BGI_T5_BGI',2916,1193,'MantaDEL')
ilm_nova <- line_dat('SV/Quartet_DNA_ILM_XTen_WUX_LCL5_4_20180703_manta_all.txt','SV/Quartet_DNA_ILM_XTen_WUX_LCL5_4_M.pos.txt','ILM_XTen',3845,973,'MantaDEL')

del_big_dat <- data.frame(rbind(bgi_500,bgi_t5,ilm_nova))
del_big_dat$type <- 'del'



bgi_500 <- line_dat('SV/Quartet_DNA_BGI_SEQ500_BGI_LCL5_1_20180328_manta_all.txt','SV/Quartet_DNA_BGI_SEQ500_BGI_LCL5_1_M.pos.txt','BGI_SEQ500_BGI',1170,626,'MantaINS')
bgi_t5 <- line_dat('SV/Quartet_DNA_BGI_T5_SZ_LCL5_1_20200716_manta_all.txt','SV/Quartet_DNA_BGI_T5_SZ_LCL5_1_M.pos.txt','BGI_T5_BGI',984,1398,'MantaINS')
ilm_nova <- line_dat('SV/Quartet_DNA_ILM_XTen_WUX_LCL5_4_20180703_manta_all.txt','SV/Quartet_DNA_ILM_XTen_WUX_LCL5_4_M.pos.txt','ILM_XTen',1530,1145,'MantaINS')

ins_big_dat <- data.frame(rbind(bgi_500,bgi_t5,ilm_nova))
ins_big_dat$type <- 'ins'


sv_big_dat <- data.frame(rbind(del_big_dat,ins_big_dat))
sv_big_dat$group <- paste(sv_big_dat$site,sv_big_dat$type,sep="_")


p <- ggplot(data=big_dat, aes(x=mendelian_violation, y=mendelian_consistent, group=group)) +
  #geom_point()+
  geom_line(aes(color=site,linetype=type))+
  #scale_color_manual(values=c('#c64b2b','#f9ce8c','#c3e3e5','#589bad','#1a1d1e','#748b42'))+
  theme_light()+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))
  coord_cartesian(ylim=c(0.9,1.0),xlim=c(0.75,1.0))


pdf('sv_filtration_threshold.pdf',height=2.5,width = 4)
plot(p)
dev.off()




###snv
all_table = 'snv_pos/BGI_SEQ500_BGI_LCL5_1_info.txt'
type = 'INDEL'

line_dat <-function(all_table,tag,snv_mc_number,snv_mv_number,indel_mc_number,indel_mv_number){
  bgi_500_all <- read.csv(all_table,header=F)
  bgi_500_all$V6[which(bgi_500_all$V6 == "")] <- 'MV'
  bgi_500_all$V6[which(bgi_500_all$V6 != "MV")] <- 'MC'
  bgi_500_all$V7 <- 'SNV'
  bgi_500_all$V7[which((nchar(bgi_500_all$V3) > 1)|(nchar(bgi_500_all$V4) > 1))] <- 'INDEL'
  ###SNV
  snv <- bgi_500_all[which(bgi_500_all$V7=='SNV'),]
  dat_point <- seq(1,2000,10)
  result_snv <- matrix(NA,ncol=2,nrow = 200)
  for(i in 1:length(dat_point)){
    print(i)
    loc <- which(snv$V5 > dat_point[i])
    sub_men <- snv[loc,]
    result_snv[i,1] <- length(which(sub_men$V6 == 'MC'))
    result_snv[i,2] <- length(which(sub_men$V6 == 'MV'))
  }
  result_snv[,1] <- result_snv[,1]/snv_mc_number
  result_snv[,2] <- result_snv[,2]/snv_mv_number
  colnames(result_snv) <- c('mendelian_consistent','mendelian_violation')
  result_snv <- data.frame(result_snv)
  result_snv$site <- tag
  result_snv$para <- dat_point
  result_snv$type <- 'SNV'

  
  indel <- bgi_500_all[which(bgi_500_all$V7=='INDEL'),]
  dat_point <- seq(1,2000,10)
  result_indel <- matrix(NA,ncol=2,nrow = 200)
  for(i in 1:length(dat_point)){
    print(i)
    loc <- which(indel$V5 > dat_point[i])
    sub_men <- indel[loc,]
    result_indel[i,1] <- length(which(sub_men$V6 == 'MC'))
    result_indel[i,2] <- length(which(sub_men$V6 == 'MV'))
  }
  result_indel[,1] <- result_indel[,1]/indel_mc_number
  result_indel[,2] <- result_indel[,2]/indel_mv_number
  colnames(result_indel) <- c('mendelian_consistent','mendelian_violation')
  result_indel <- data.frame(result_indel)
  result_indel$site <- tag
  result_indel$para <- dat_point 
  result_indel$type <- 'INDEL'
  
  result <- list(result_snv,result_indel)
  return(result)
}

bgi_500 <- line_dat('snv_pos/BGI_SEQ500_BGI_LCL5_1_info.txt','BGI_SEQ500_BGI',3781643,111662,647683,237430)
bgi_t5 <- line_dat('snv_pos/BGI_SEQT5_BSZ_LCL5_1_info.txt','BGI_T5_BGI',3805528,118565,709030,218940)
ilm_xten <- line_dat('snv_pos/ILM_Nova_WUX_LCL5_4_info.txt','ILM_XTen',3821217,135036,676516,208870)


snv_dat <- data.frame(rbind(bgi_500[[1]],bgi_t5[[1]],ilm_xten[[1]]))
indel_dat <- data.frame(rbind(bgi_500[[2]],bgi_t5[[2]],ilm_xten[[2]]))
small_variants_big_dat <- data.frame(rbind(snv_dat,indel_dat))
small_variants_big_dat$group <- paste(small_variants_big_dat$site,small_variants_big_dat$type,sep="_")

big_dat <- data.frame(rbind(sv_big_dat,small_variants_big_dat))
p <- ggplot(data=big_dat, aes(x=mendelian_violation, y=mendelian_consistent, group=group)) +
  #geom_point()+
  geom_line(aes(color=site,linetype=type))+
  scale_linetype_manual(values=c("solid","dotted",'dashed','twodash'))+
  #scale_color_manual(values=c('#c64b2b','#f9ce8c','#c3e3e5','#589bad','#1a1d1e','#748b42'))+
  theme_light()+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))
coord_cartesian(ylim=c(0.9,1.0),xlim=c(0.75,1.0))


pdf('filtration_threshold.pdf',height=2.5,width = 4.5)
plot(p)
dev.off()
