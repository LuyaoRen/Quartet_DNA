ibrary(reshape2)
indel <- read.table('small_variants_indel_jaccard_index.txt',header = T)
indel_sub <- indel[- grep('RPG',rownames(indel)),-grep('RPG',colnames(indel))]
indel_sub <- indel[grep('LCL',rownames(indel)),grep('LCL',colnames(indel))]
indel_sub$sample1 <- rownames(indel_sub)
indel_sub$library <- sapply(strsplit(indel_sub$sample1,"_"),function(x){x[9]})
long_indel <- melt(indel_sub)

tech1 <- sapply(strsplit(long_indel$sample1,"_"),function(x){x[3]})
platform1 <- sapply(strsplit(long_indel$sample1,"_"),function(x){x[4]})
site1 <- sapply(strsplit(long_indel$sample1,"_"),function(x){x[5]})
sample1 <- sapply(strsplit(long_indel$sample1,"_"),function(x){x[6]})
date1 <- sapply(strsplit(long_indel$sample1,"_"),function(x){x[8]})

idenifier1 <- paste(tech1,platform1,site1,sample1,date1,sep="_")


long_indel$variable <- as.character(long_indel$variable)
tech2 <- sapply(strsplit(long_indel$variable,"_"),function(x){x[3]})
platform2 <- sapply(strsplit(long_indel$variable,"_"),function(x){x[4]})
site2 <- sapply(strsplit(long_indel$variable,"_"),function(x){x[5]})
sample2 <- sapply(strsplit(long_indel$variable,"_"),function(x){x[6]})
date2 <- sapply(strsplit(long_indel$variable,"_"),function(x){x[8]})

idenifier2 <- paste(tech2,platform2,site2,sample2,date2,sep="_")

indel_rep <- long_indel[idenifier1==idenifier2,]
indel_rep <- indel_rep[which(indel_rep$value != 1),]

####
long_indel$type <- "Complex"
long_indel$type[which((tech1 != tech2)&(sample1 == sample2))] <- "Tech"
long_indel$type[which((tech1 == tech2)&(platform1 != platform2)&(site1 == site2)&(sample1 == sample2))] <- "Platform"
long_indel$type[which((tech1 == tech2)&(platform1 == platform2)&(site1 != site2)&(sample1 == sample2))] <- "Site"
long_indel$type[which(idenifier1==idenifier2)] <- "Intra-Btach"



tech <- sapply(strsplit(indel_rep$sample1,"_"),function(x){x[3]})
platform <- sapply(strsplit(indel_rep$sample1,"_"),function(x){x[4]})
site <- sapply(strsplit(indel_rep$sample1,"_"),function(x){x[5]})
date <- sapply(strsplit(indel_rep$sample1,"_"),function(x){x[8]})

indel_rep$batch <- paste(tech,platform,site,date,sep="_")


snv <- read.table('small_variants_snv_jaccard_index.txt',header = T)
snv_sub <- snv[- grep('RPG',rownames(snv)),-grep('RPG',colnames(snv))]
snv_sub <- snv[grep('LCL',rownames(snv)),grep('LCL',colnames(snv))]
snv_sub$sample1 <- rownames(snv_sub)
snv_sub$library <- sapply(strsplit(snv_sub$sample1,"_"),function(x){x[9]})
long_snv <- melt(snv_sub)

tech1 <- sapply(strsplit(long_snv$sample1,"_"),function(x){x[3]})
platform1 <- sapply(strsplit(long_snv$sample1,"_"),function(x){x[4]})
site1 <- sapply(strsplit(long_snv$sample1,"_"),function(x){x[5]})
sample1 <- sapply(strsplit(long_snv$sample1,"_"),function(x){x[6]})
date1 <- sapply(strsplit(long_snv$sample1,"_"),function(x){x[8]})

idenifier1 <- paste(tech1,platform1,site1,sample1,date1,sep="_")


long_snv$variable <- as.character(long_snv$variable)
tech2 <- sapply(strsplit(long_snv$variable,"_"),function(x){x[3]})
platform2 <- sapply(strsplit(long_snv$variable,"_"),function(x){x[4]})
site2 <- sapply(strsplit(long_snv$variable,"_"),function(x){x[5]})
sample2 <- sapply(strsplit(long_snv$variable,"_"),function(x){x[6]})
date2 <- sapply(strsplit(long_snv$variable,"_"),function(x){x[8]})

idenifier2 <- paste(tech2,platform2,site2,sample2,date2,sep="_")

snv_rep <- long_snv[idenifier1==idenifier2,]
snv_rep <- snv_rep[which(snv_rep$value != 1),]

tech <- sapply(strsplit(snv_rep$sample1,"_"),function(x){x[3]})
platform <- sapply(strsplit(snv_rep$sample1,"_"),function(x){x[4]})
site <- sapply(strsplit(snv_rep$sample1,"_"),function(x){x[5]})
date <- sapply(strsplit(snv_rep$sample1,"_"),function(x){x[8]})

snv_rep$batch <- paste(tech,platform,site,date,sep="_")

indel_rep$type <- "INDEL"
snv_rep$type <- "SNV"

dat <- data.frame(rbind(indel_rep,snv_rep))

p <- ggplot(dat, aes(x=batch, y=value,fill=library)) + 
  geom_boxplot()+
  ylab(c('Jaccard Index'))+
  xlab(c(''))+
  facet_grid(type ~ library,scales="free")+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust = 1 ,
                                   size=12),
        axis.text.y = element_text( size=12))+
  theme(
    axis.title.x = element_text( size=14),
    axis.title.y = element_text( size=14)
  )
pdf('Jaccard_Index.pdf',height = 6.5,width  =7)
plot(p)
dev.off()

#############################
# 不同数据集间检测突变的一致性
#############################
library(ComplexHeatmap)
library(circlize)
snv <- read.table('small_variants_snv_jaccard_index.txt')
indel <- read.table('small_variants_indel_jaccard_index.txt')

snv_lcl5 <- snv[grep('LCL8',rownames(snv)),grep('LCL8',colnames(snv))]
indel_lcl5 <- indel[grep('LCL8',rownames(indel)),grep('LCL8',colnames(indel))]

snv_lcl5_sorted <- snv_lcl5[ order(colnames(snv_lcl5)), order(colnames(snv_lcl5))]
indel_lcl5_sorted <- indel_lcl5[ order(colnames(indel_lcl5)), order(colnames(indel_lcl5))]

lcl5_pd <- data.frame(colnames(indel_lcl5_sorted))
lcl5_pd$platform <- sapply(strsplit(lcl5_pd[,1],"_"),function(x){x[4]})
lcl5_pd$site <- sapply(strsplit(lcl5_pd[,1],"_"),function(x){x[5]})
lcl5_pd$library <- c('PCR','PCR','PCR','PCR-free','PCR-free','PCR-free',
                     'PCR-free','PCR-free','PCR-free','PCR-free','PCR-free','PCR-free',
                     'PCR-free','PCR-free','PCR-free','PCR','PCR','PCR',
                     'PCR','PCR','PCR','PCR','PCR','PCR','PCR','PCR','PCR')
colnames(lcl5_pd)[1] <- 'sample'
lcl5_pd$coverage <- c(48.15,46.72,48.03,
                      36.3,33.71,36.63,
                      40.21,31.65,42.7,52.52,56.64,58.99,
                      63.42,47.34,48.97,
                      36.43,33.79,29.27,
                      26.49,26.38,25.65,
                      28.33,28.64,28.51,
                      39.31,41.49,36.29)



ji_col_fun = colorRamp2(c(0.5, 0.75, 1), c("#377EB8", "white", "#E41A1C"))
cov_col_fun = colorRamp2(c(20, 45, 70), c("green", "yellow", "red"))
ha = HeatmapAnnotation(df=lcl5_pd[,2:5],
                       col = list(platform = c("Nova" = "#7570B3", "SEQ2000" = "#E7298A","T7"="#1B9E77","XTen" = "#D95F02"),
                                  site = c("ARD"="#4DAF4A","BGI"="#FFFF33",'BRG'='#377EB8','NVG'='#984EA3','WGE'='#E41A1C','WUX'='#FF7F00'),
                                  library = c('PCR' = '#8e0c24','PCR-free' = '#587aa5'),
                                  coverage = cov_col_fun))
a <- Heatmap(as.matrix(indel_lcl5_sorted),top_annotation = ha, col = ji_col_fun, show_row_names = FALSE,show_column_names = FALSE, 
             show_column_dend = FALSE, column_title = "INDEL",name = "Jaccard Index")

pdf('LCL8_indel_ji.pdf',height = 4,width = 5)
a
dev.off()



snv_lcl <- snv[grep('LCL',rownames(snv)),grep('LCL',colnames(snv))]
indel_lcl <- indel[grep('LCL',rownames(indel)),grep('LCL',colnames(indel))]

snv_lcl_sorted <- snv_lcl5[ order(colnames(snv_lcl)), order(colnames(snv_lcl))]
indel_lcl_sorted <- indel_lcl5[ order(colnames(indel_lcl)), order(colnames(indel_lcl))]

lcl_pd <- data.frame(rownames(indel_lcl_sorted))
lcl_pd$platform <- sapply(strsplit(lcl_pd[,1],"_"),function(x){x[4]})
lcl_pd$site <- sapply(strsplit(lcl_pd[,1],"_"),function(x){x[5]})
lcl_pd$library <- sapply(strsplit(lcl_pd[,1],"_"),function(x){x[9]})
lcl_pd$sample <- sapply(strsplit(lcl_pd[,1],"_"),function(x){x[6]})
colnames(lcl_pd)[1] <- 'id'


ji_col_fun = colorRamp2(c(0.2, 0.6, 1), c("#377EB8", "white", "#E41A1C"))
ha = HeatmapAnnotation(df=lcl_pd[,2:5],
                       col = list(platform = c("Nova" = "#7570B3", "SEQ2000" = "#E7298A","T7"="#1B9E77","XTen" = "#D95F02"),
                                  site = c("ARD"="#4DAF4A","BGI"="#FFFF33",'BRG'='#377EB8','NVG'='#984EA3','WGE'='#E41A1C','WUX'='#FF7F00'),
                                  library = c('PCR' = '#8e0c24','pcr-free' = '#587aa5'),
                                  sample = c("LCL8" = "#F16745","LCL7"="#FFC65D","LCL6"="#7BC8A4","LCL5"="#4CC3D9")))
a <- Heatmap(as.matrix(indel_lcl_sorted),top_annotation = ha, col = ji_col_fun, show_row_names = FALSE,show_column_names = FALSE, 
             show_column_dend = FALSE, column_title = "INDEL",name = "Jaccard Index")

pdf('INDEL_ji.pdf',height = 4,width = 5.5)
a
dev.off()

