
rm(list = ls())

library(MetaboAnalystR)
library(ropls)
library(ggplot2)
library(dplyr)
library(plyr)
library(openxlsx)

outdir <- "D:/R_space/Work_space/pipeline-5"
subdir <- c("1.Information","2.QualityControl","3.Statistics","4.SignificanceAnalysis","5.PathwayAnalysis")

for (i in 1:length(subdir)){
  dir.create(file.path(outdir, subdir[i]))
}
setwd(outdir)
infile <- "1.Information/neg-raw-metaboAnalystInput-2.csv"
sample_infor <- read.xlsx(paste(outdir,"1.Information/Sample_Information.xlsx",sep = "/"),sheet = 1)

# if (file.exists(subdir)){
#   setwd(file.path(outdir, subdir))
# } else {
#   dir.create(file.path(outdir, subdir))
#   setwd(file.path(outdir, subdir))
# }

rawdata <- read.csv(infile,header = T,sep = ",",quote = " ",stringsAsFactors = F)
head(rawdata)
rawdata2 <- as.data.frame(lapply(rawdata[-1,-1],as.numeric))
rownames(rawdata2) <- rawdata[-1,1]
# Obtain QC sample data matrix
qcdata <- rawdata2[,grepl("QC|qc",colnames(rawdata2))] 
# Calculate RSDs and data filtering
qcmean <- apply(qcdata, 1, mean)
qcsd <- apply(qcdata, 1, sd)
qccv <- data.frame(RSD = qcsd/qcmean)
qccv$MZ <- rownames(qccv)
rsdfilter <- subset(qccv,RSD < 0.3)
head(rsdfilter)
rsdfilter <- merge(rawdata,rsdfilter,by.x = colnames(rawdata)[1],by.y = "MZ",sort = F)
rsdfilter <- rsdfilter[,!grepl("RSD",colnames(rsdfilter))]
RawFilter <- rbind(rawdata[1,],rsdfilter)
# Save data filtered
write.table(RawFilter,"data_filter_rsd.csv",sep = ",",row.names = F)


mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "./data_filter_rsd.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)
# Data Processing
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeMissingVar(mSet, method="min")
mSet<-SanityCheckData(mSet)
mSet<-FilterVariable(mSet, "iqr", "F", 25)  #### 
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-SaveTransformedData(mSet)

normalize_data <- read.csv(paste(outdir, "data_normalized.csv",sep="/"),sep = ",",header = T)

colnames(normalize_data)[1] <- "Compound ID"  
rownames(normalize_data) <- normalize_data[,1]
nor_data <- t(normalize_data)
nor_data<-data.frame(nor_data,stringsAsFactors = F) 
colnames(nor_data) <- nor_data[1,]
nor_data <- nor_data[-1,-1]
####
nor_data1 <- as.data.frame(lapply(nor_data,as.numeric)) # 数据转换成数值形式
rownames(nor_data1) <- rownames(nor_data)
colnames(nor_data1) <- colnames(nor_data)
write.table(nor_data1,"data_nor_clean.csv",sep=",", row.names = F,col.names = T)
########################################################### 
# QC PCA analysis(ropls)
########################################################### 
nor.pca <- opls(x = nor_data1,predI=4)
names(attributes(nor.pca))    #查看包含信息
head(nor.pca@scoreMN)    #例如，提取对象得分（各样本在 PCA 轴中的排序坐标）
pca.scoreMN <- as.data.frame(nor.pca@scoreMN)
pca.scoreMN$Sample <- rownames(pca.scoreMN)
# Data of PCA plot
pca.scoreMN <- merge(pca.scoreMN, sample_infor[,3:4],by.x = "Sample",by.y = "Sample.analysis.name")
head(pca.scoreMN)
write.csv(pca.scoreMN,"pca_plot_data.csv",row.names = F)
########################################################### 
# ggplot2 QC PCA analysis
########################################################### 
xlab <- paste0( "PC1(",round(nor.pca@modelDF$R2X[1]* 100, 2), "%)")
ylab <- paste0( "PC2(",round(nor.pca@modelDF$R2X[2]* 100, 2), "%)")
p1 <- ggplot(data = pca.scoreMN,aes(x=p1,y=p2,color = Group ))+
  stat_ellipse(aes(fill=Group),type = "norm", geom = "polygon",alpha= 0.3,color=NA)+
  geom_point(size=3) + labs(x=xlab,y=ylab,size=5) + guides(fill="none")+
  theme_bw() +theme(panel.grid=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  geom_hline(yintercept = 0,size=0.7)+geom_vline(xintercept = 0,size=0.7)+
  theme(legend.title = element_blank())+
  theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+
  theme(axis.title.y = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  ggtitle("QC-PCA")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) #设置标题居中

p1

ggsave(p1,filename = "2.QualityControl/QC-PCA.pdf",width = 12,height = 9,dpi = 300)
ggsave(p1,filename = "2.QualityControl/QC-PCA.png",width = 12,height = 9,dpi = 300)

##########################################################
# QC heatmap
##########################################################
library(pheatmap)
library(reshape2)

data_heat <- RawFilter[-1,-1]
rownames(sample_infor) <- sample_infor$Sample.analysis.name
Group <-data.frame(class= sample_infor[,4])
rownames(Group) <- rownames(sample_infor)
df_heat <- as.data.frame(lapply(data_heat,as.numeric))

p2 <- pheatmap(df_heat,scale = "row", cluster_cols = T,show_rownames = F,cluster_rows = T,annotation_col = Group,legend_labels = NA,
              show_colnames=T,col=colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 30, treeheight_col = 30)
p2
ggsave(p2,filename = "2.QualityControl/QC-heatmap.pdf",width = 12,height = 9,dpi = 300)
ggsave(p2,filename = "2.QualityControl/QC-heatmap.png",width = 12,height = 9,dpi = 300)

# QC boxplot
##########################################################
library(tidyr)

databox <- log10(df_heat)
databox2 <- gather(databox,key = "Sample",value = "value")
databox2 <- merge(databox2,sample_infor[,3:4],by.x ="Sample",by.y = "Sample.analysis.name",all.x = T)

b <- ggplot(databox2,aes(factor(Sample),value, fill=Group)) + geom_boxplot()+ 
  theme(axis.text = element_text(size=7,angle = 90), axis.title = element_text(size=12)) +
  labs(color="Group", x="Sample",y="log10(intensity)")
b
ggsave(b,filename = "2.QualityControl/QC-IntensityBoxplot.pdf",width = 12,height = 9,dpi = 300)
ggsave(b,filename = "2.QualityControl/QC-IntensityBoxplot.png",width = 12,height = 9,dpi = 300)

################# QualityControl END ####################

################ Statistics Start #######################

mSet<-PlotNormSummary(mSet, "3.Statistics/norm_0_", "png", 300, width=NA)
mSet<-PlotSampleNormSummary(mSet, "3.Statistics/snorm_0_", "png", 300, width=NA)
mSet<-SaveTransformedData(mSet)  ## Save processed data

#################################################################
#Step 1:Statistic analysis
#################################################################
mSet<-FC.Anal(mSet, 2.0, 1, FALSE) # 1 或者 0调整比较组的前后顺序
mSet<-PlotFC(mSet, "3.Statistics/fc_2_", "png", 300, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "raw", FALSE)
mSet<-PlotTT(mSet, "3.Statistics/tt_raw_", "png", 300, width=NA)

# mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
# mSet<-PlotTT(mSet, "tt_fdr_", "png", 300, width=NA)

Fc_tt <- as.data.frame(cbind(mSet$analSet$tt$p.value,mSet$analSet$fc$fc.all))
colnames(Fc_tt) <- c("P.value","FC")
Fc_tt$`Compound ID` <- rownames(Fc_tt)

cut_off_pvalue = 0.05
cut_off_FC = 2
Fc_tt$Change[(Fc_tt$P.value > cut_off_pvalue)|(1/cut_off_FC <= Fc_tt$FC)& Fc_tt$FC <= cut_off_FC] <- "None"
Fc_tt$Change[(Fc_tt$P.value < cut_off_pvalue)&(cut_off_FC < Fc_tt$FC)] <- "Up"
Fc_tt$Change[(Fc_tt$P.value < cut_off_pvalue)&(Fc_tt$FC < 1/cut_off_FC)] <- "Down"
aa <- Fc_tt %>% select('Compound ID','P.value','FC','Change',everything())
write.csv(aa,"3.Statistics/FoldChange-Pvalue.csv",row.names = F)
# 绘制火山图====================================
x1 <- log2(Fc_tt$FC)
y1 <- -log10(Fc_tt$P.value)

p <- ggplot(
  #设置数据
  Fc_tt, 
  aes(x = x1, 
      y = y1, 
      colour=Change)) +
  geom_point(alpha=0.4, size=2.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  xlim(-5,5)+
  # 坐标轴
  labs(x="log2(Fold Change)",
       y="-log10 (P-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )
p
ggsave(p,filename = "3.Statistics/Volcano.pdf",width = 12,height = 9,dpi = 300)
ggsave(p,filename = "3.Statistics/Volcano.png",width = 12,height = 9,dpi = 300)


# ropls packages statistic

########################################################### 
# PCA analysis(ropls)
###########################################################
nor_data2 <- nor_data1[!grepl("QC|qc|Qc|qC",rownames(nor_data1)),]

nor.pca <- opls(x = nor_data2)
names(attributes(nor.pca))    #查看包含信息
head(nor.pca@scoreMN)    #例如，提取对象得分（各样本在 PCA 轴中的排序坐标）
# pca.scoreMN <- nor.pca@scoreMN
# pca.scoreMN <- cbind(pca.scoreMN, nor_sample)

pca.scoreMN <- as.data.frame(nor.pca@scoreMN)
pca.scoreMN$Sample <- rownames(pca.scoreMN)
pca.scoreMN <- merge(pca.scoreMN, sample_infor[,3:4],by.x = "Sample",by.y = "Sample.analysis.name")
head(pca.scoreMN)

########################################################### 
# ggplot2 PCA analysis
########################################################### 
xlab <- paste0( "PC1(",round(nor.pca@modelDF$R2X[1]* 100, 2), "%)")
ylab <- paste0( "PC2(",round(nor.pca@modelDF$R2X[2]* 100, 2), "%)")
p1 <- ggplot(data = pca.scoreMN,aes(x=p1,y=p2,color = Group ))+
  stat_ellipse(aes(fill=Group),type = "norm", geom = "polygon",alpha= 0.3,color=NA)+
  geom_point(size=3) + labs(x=xlab,y=ylab,size=5) + guides(fill="none")+
  theme_bw() +theme(panel.grid=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  geom_hline(yintercept = 0,size=0.7)+geom_vline(xintercept = 0,size=0.7)+
  theme(legend.title = element_blank())+
  theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+
  theme(axis.title.y = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  ggtitle("PCA")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) #设置标题居中

p1

ggsave(p1,filename = "3.Statistics/PCA.pdf",width = 12,height = 9,dpi = 300)
ggsave(p1,filename = "3.Statistics/PCA.png",width = 12,height = 9,dpi = 300)

###########################################################
## PLS-DA analysis
###########################################################
nor.plsda <- opls(x = nor_data2, y = sample_infor[!grepl("QC|qc|Qc|qC",sample_infor$Group), 'Group'], orthoI = 0)
names(attributes(nor.plsda))
head(nor.plsda@scoreMN)    #例如，样本在 PLS-DA 轴上的位置
head(nor.plsda@loadingMN)    #例如，代谢物在 PLS-DA 轴上的位置

plsda.scoreMN <- as.data.frame(nor.plsda@scoreMN)
plsda.scoreMN$Sample <- rownames(plsda.scoreMN)
plsda.scoreMN <- merge(plsda.scoreMN, sample_infor[,3:4],by.x = "Sample",by.y = "Sample.analysis.name")
head(plsda.scoreMN)
#plsda.scoreMN <- cbind(plsda.scoreMN, nor_sample)


#########################################################
# ggplot2 PLS-DA analysis
###########################################################
xlab1 <- paste0( "PC1(",round(nor.plsda@modelDF$R2X[1]* 100, 2), "%)")
ylab1 <- paste0( "PC2(",round(nor.plsda@modelDF$R2X[2]* 100, 2), "%)")
p1 <- ggplot(data = plsda.scoreMN,aes(x=p1,y=p2,color = Group ))+
  stat_ellipse(aes(fill=Group),type = "norm", geom = "polygon",alpha= 0.3,color=NA)+
  geom_point(size=3) + labs(x=xlab1,y=ylab1,size=5) + guides(fill="none")+
  theme_bw() +theme(panel.grid=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  geom_hline(yintercept = 0,size=0.7)+geom_vline(xintercept = 0,size=0.7)+
  theme(legend.title = element_blank())+
  theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+
  theme(axis.title.y = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  ggtitle("PLS-DA") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) #设置标题居中

p1
ggsave(p1,filename = "3.Statistics/PLS-DA.pdf",width = 12,height = 9,dpi = 300)
ggsave(p1,filename = "3.Statistics/PLS-DA.png",width = 12,height = 9,dpi = 300)

#plsda.permu <- nor.plsda@suppLs$permMN # plsda permutation data

###########################################################
# OPLS-DA analysis
###########################################################
nor.oplsda <- opls(x = nor_data2, y = sample_infor[!grepl("QC|qc|Qc|qC",sample_infor$Group), 'Group'],predI=1,orthoI=NA,permI=200) #orthoI = NA 时执行 OPLS-DA
names(attributes(nor.oplsda))
head(nor.oplsda@scoreMN)    #例如，样本在 OPLS-DA 轴上的位置
head(nor.oplsda@loadingMN)    #例如，代谢物在 OPLS-DA 轴上的位置

oplsda.scoreMN <- nor.oplsda@scoreMN
oplsda.orthoScoreMN <- nor.oplsda@orthoScoreMN
oplsda.score <- as.data.frame(cbind(nor.oplsda@scoreMN,nor.oplsda@orthoScoreMN))
oplsda.score$Sample <- rownames(oplsda.score)
oplsda.score <- merge(oplsda.score,sample_infor[,3:4],by.x = "Sample",by.y = "Sample.analysis.name")
#oplsda.score <- cbind(nor.oplsda@scoreMN,nor.oplsda@orthoScoreMN, nor_sample)


# ggplot2 OPLS-DA analysis
xlab1 <- paste0( "PC1(",round(nor.oplsda@modelDF$R2X[1]* 100, 2), "%)")
ylab1 <- paste0( "PCo1(",round(nor.oplsda@modelDF$R2X[2]* 100, 2), "%)")
p1 <- ggplot(data = oplsda.score,aes(x=p1,y=o1,color = Group ))+
  stat_ellipse(aes(fill=Group),type = "norm", geom = "polygon",alpha= 0.3,color=NA)+
  geom_point(size=3) + labs(x=xlab1,y=ylab1,size=5) + guides(fill="none")+
  theme_bw() +theme(panel.grid=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  geom_hline(yintercept = 0,size=0.7)+geom_vline(xintercept = 0,size=0.7)+
  theme(legend.title = element_blank())+
  theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+
  theme(axis.title.y = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  ggtitle("OPLS-DA") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) #设置标题居中

p1
ggsave(p1,filename = "3.Statistics/OPLS-DA.pdf",width = 12,height = 9,dpi = 300)
ggsave(p1,filename = "3.Statistics/OPLS-DA.png",width = 12,height = 9,dpi = 300)

###########################################################
# OPLS-DA loading analysis
###########################################################
# opls.loading <- cbind(nor.oplsda@loadingMN,nor.oplsda@orthoLoadingMN)
# opls.loading <- as.data.frame(opls.loading)
# xlab1 <- paste0( "PC1(",round(nor.oplsda@modelDF$R2X[1]* 100, 2), "%)")
# ylab1 <- paste0( "PCo1(",round(nor.oplsda@modelDF$R2X[2]* 100, 2), "%)")
# p1 <- ggplot(data = opls.loading,aes(x=p1,y=o1 ))+
#     geom_point(size=1.5,color="red") + labs(x=xlab1,y=ylab1,size=5) + guides(fill="none")+
#   geom_hline(yintercept = 0,size=0.7)+geom_vline(xintercept = 0,size=0.7)+
#   theme_bw() +theme(panel.grid=element_blank())+
#   theme(panel.border = element_rect(fill=NA,color="black", size=1.1, linetype="solid"))+ # 边框线
#   theme(legend.title = element_blank())+
#   theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+ # 横轴标签
#   theme(axis.title.y = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+ # 纵轴标签
#   theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+ # 横轴刻度
#   theme(axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+ # 纵轴刻度
#   ggtitle("Loading-OPLS-DA") +
#   theme(plot.title = element_text(hjust = 0.5,face = "bold")) #设置标题居中
# p1
# 
# ggsave(p1,filename = "Loading-OPLS-DA.pdf",width = 12,height = 9,dpi = 300)
# ggsave(p1,filename = "Loading-OPLS-DA.png",width = 12,height = 9,dpi = 300)
# 
# ###########################################################
# # Permutation Test
# ###########################################################
# opls.permutation <- as.data.frame(nor.oplsda@suppLs$permMN)
# png(filename = "Permutation-OPLS-DA.png")
# permuPlot <- plot(nor.oplsda, typeVc = 'permutation')
# dev.off()
# 
# write.table(opls.permutation,"OPLS-permutation.csv",sep = ",",row.names = F)
# 
# 
# 
# #VIP 值帮助寻找重要的代谢物
# vipVn <-as.data.frame(getVipVn(nor.oplsda))
# names(vipVn)[1] <- "VIP"
# aa2 <- cbind(aa,vipVn)
# aa2 <- aa2 %>% dplyr::select('Compound ID','P.value','FC','VIP','Change')
# df <- cbind(normalize_data[-1,],aa2[,-1])
# write.csv(df,"no-significant.csv",row.names = F)
# 
# 
# ###########################################################
# #Step 2 Function analysis
# ##########################################################
# outdir3 <- "D:/R_space/Work_space/pipeline-3/Function-analysis/"
# setwd(outdir3)
# data_Funct <- InitDataObjects("mass_table", "mummichog", FALSE)
# data_Funct<-SetPeakFormat(data_Funct, "colu")
# data_Funct<-UpdateInstrumentParameters(data_Funct, 5.0, "negative", "yes", 0.02);
# data_Funct<-SetRTincluded(data_Funct, "no")
# data_Funct<-Read.TextData(data_Funct, "../FC-Ttest-Analysis/data_processed.csv", "mpt", "disc");
# data_Funct<-SanityCheckMummichogData(data_Funct)
# data_Funct<-ReplaceMin(data_Funct);
# data_Funct<-SanityCheckMummichogData(data_Funct)
# data_Funct<-FilterVariable(data_Funct, "none", "F", 25)
# data_Funct<-PreparePrenormData(data_Funct)
# data_Funct<-Normalization(data_Funct, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
# data_Funct<-PlotNormSummary(data_Funct, "norm_function_", "png", 72, width=NA)
# data_Funct<-PlotSampleNormSummary(data_Funct, "snorm_function_", "png", 72, width=NA)
# add.vec <- c("M-H [1-]","M-2H [2-]","M-H2O-H [1-]","M-H+O [1-]","M+K-2H [1-]","M+Na-2H [1- ]","M+Cl [1-]","M+CH3COO [1-]")
# data_Funct<-Setup.AdductData(data_Funct, add.vec);
# data_Funct<-PerformAdductMapping(data_Funct, "negative")
# data_Funct<-SetPeakEnrichMethod(data_Funct, "mum", "v2")
# ###########################################################
# # data_Funct<-DoPeakConversion(data_Funct)         # ?  DoPeakConversion 功能缺失          
# # data_Funct<-SetMummichogPval(data_Funct, 1.0E-4) # ? P-value cutoff
# # data_Funct<-PerformPSEA(data_Funct, "hsa_kegg", "current", 3 , 100)
# # data_Funct<-PlotPeaks2Paths(data_Funct, "peaks_to_paths_0_", "png", 72, width=NA)
# 
# 
# ######################################################
# #Step 3 Pathway analysis
# ######################################################
# outdir4 <- "D:/R_space/Work_space/pipeline-3/Pathway-Analysis/"
# setwd(outdir4)
# infile2 <- read.delim(paste(outdir3,"mummichog_matched_compound_all.csv",sep="/"),header = T,sep = ",")
# keggID <- infile2$Matched.Compoun
# mSet<-InitDataObjects("conc", "pathora", FALSE)
# mSet<-Setup.MapData(mSet, keggID[1:20]);
# mSet<-CrossReferencing(mSet, "kegg");
# mSet<-CreateMappingResultTable(mSet)
# mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
# mSet<-SetMetabolomeFilter(mSet, F);
# mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
# mSet<-PlotPathSummary(mSet, T, "path_view_0_", "png", 72, width=NA)
# mSet<-SaveTransformedData(mSet)
# 
# #######################################################
# #Step 4 Enrichment analysis
# #######################################################
# setwd("D:/R_space/Work_space/pipeline-3/Pathway-Enrichment/")
# mSet<-InitDataObjects("conc", "msetora", FALSE)
# mSet<-Setup.MapData(mSet, keggID[1:20]);
# mSet<-CrossReferencing(mSet, "kegg");
# mSet<-CreateMappingResultTable(mSet)
# mSet<-SetMetabolomeFilter(mSet, F);
# mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
# mSet<-CalculateHyperScore(mSet)
# mSet<-PlotORA(mSet, "ora_0_", "net", "png", 72, width=NA)
# mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 72, width=NA)



