
#=====================================================================================
#
#  Code chunk 1 CGGA
#
#=====================================================================================


##### CGGA data
Survival_2=read.delim("CGGA.mRNAseq_693.clinical.20190701.txt")
Survival_2=Survival_2[,c(1,8,7,6,2)]
Survival_2=Survival_2[Survival_2$PRS_type=="Primary",]
Survival_2=na.omit(Survival_2)
rownames(Survival_2)=Survival_2[,1]
Survival_2=Survival_2[,c(1:4)]
colnames(Survival_2)=c('SampleID','OS_event','OS_time','Age')

Gene_2=read.delim("CGGA.mRNAseq_693.RSEM-genes.20190701.txt")
Gene_2=as.matrix(Gene_2)
rownames(Gene_2)=Gene_2[,1]
Gene_2=Gene_2[,-1]  # 23987   693

Gene_2=Gene_2[,intersect(colnames(Gene_2),rownames(Survival_2))]            
Survival_2=Survival_2[intersect(colnames(Gene_2),rownames(Survival_2)),] 


#=====================================================================================
#
#  Code chunk 2 TCGA
#
#=====================================================================================


##### TCGA GBM
Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,56,106)]
Survival_GBM=Survival_GBM[Survival_GBM$sample_type=="Primary Tumor",]
Survival_GBM=na.omit(Survival_GBM)   
rownames(Survival_GBM)=Survival_GBM[,1]
Survival_GBM=Survival_GBM[,c(1:4)]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','Age')

Gene_GBM=read.csv("TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Gene_GBM=Gene_GBM[,intersect(colnames(Gene_GBM),rownames(Survival_GBM))]
Survival_GBM=Survival_GBM[intersect(colnames(Gene_GBM),rownames(Survival_GBM)),]


##### TCGA LGG
Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,20,87)]
Survival_LGG=Survival_LGG[Survival_LGG$sample_type=="Primary Tumor",]
Survival_LGG=na.omit(Survival_LGG)   
rownames(Survival_LGG)=Survival_LGG[,1]
Survival_LGG=Survival_LGG[,c(1:4)]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','Age')

Gene_LGG=read.csv("TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,1]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-1]

Gene_LGG=Gene_LGG[,intersect(colnames(Gene_LGG),rownames(Survival_LGG))]
Survival_LGG=Survival_LGG[intersect(colnames(Gene_LGG),rownames(Survival_LGG)),]

##### TCGA
Gene_TCGA=cbind(Gene_GBM[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),],Gene_LGG[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),])

Survival_TCGA=rbind(Survival_GBM,Survival_LGG)


#=====================================================================================
#
#  Code chunk 3 GEO
#
#=====================================================================================


Survival_geo=read.csv("GSE55918-clinical information.csv")
Survival_geo=Survival_geo[,c(1,14,13)]                      ##n.sample=1841
rownames(Survival_geo)=Survival_geo[,1]
colnames(Survival_geo)=c('SampleID','OS_event','OS_time')
Survival_geo=Survival_geo[Survival_geo$OS_time != "Null", ] #n.sample=1408

Gene_geo=read.csv("GSE55918_GeneExpressionData_Microarray.csv")
Gene_geo=as.matrix(Gene_geo)
rownames(Gene_geo)=Gene_geo[,1]
Gene_geo=Gene_geo[,-1]  

Gene_geo=Gene_geo[,intersect(colnames(Gene_geo),rownames(Survival_geo))]            
Survival_geo=Survival_geo[intersect(colnames(Gene_geo),rownames(Survival_geo)),] #n.sample=1003


gse <- getGEO("GSE16011", GSEMatrix = TRUE)

gpl8542 = gse$GSE16011_series_matrix.txt.gz



f <- fData(gpl8542)

f[f$ID=="3965_at",]
#ID                                          Description
#3965_at 3965_at lectin, galactoside-binding, soluble, 9 (galectin 9)

f[f$Description=="hepatitis A virus cellular receptor 2",]  
#84868_at 84868_at hepatitis A virus cellular receptor 2

p=pData(gpl8542)   
sample_8542=rownames(p)
sample_55918=rownames(Survival_geo)
sample_55918=str_sub(sample_55918, 1, 9)
samples=intersect(sample_8542, sample_55918)

Survival_16011=Survival_geo[samples,]

x=exprs(gpl8542)
x=x[,samples]

Survival_Gene=x[c("3965_at","84868_at"),]  #Gal9/Tim3

time=as.numeric(Survival_16011[,3])
status=as.numeric(Survival_16011[,2])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))

Gene_marker=rbind(Gene_marker,Gene_marker[1,]*Gene_marker[2,])


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


Marker_genes=c('LGALS9', 'HAVCR2')


##### data (1) CGGA
Survival_Gene=Gene_2[t(Marker_genes),]
Survival=as.matrix(Survival_2)

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))

Gene_marker=rbind(Gene_marker,Gene_marker[1,]*Gene_marker[2,],Age)

#Gene_marker=Gene_marker[1,]


##### data (2) TCGA
Survival_Gene=Gene_TCGA[t(Marker_genes),]
Survival=as.matrix(Survival_TCGA)

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))

Gene_marker=rbind(Gene_marker,Gene_marker[1,]*Gene_marker[2,],Age)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


predicted=Gene_marker[1,]

predicted0=as.numeric(predicted)
predicted0=predicted[!is.na(status)]
status0=na.omit(status)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc1=roc(status0,as.numeric(predicted0))
AUC=roc1$auc
AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc1$sensitivities,roc1$specificities)))
## optimal cut-off point 
sort(predicted0,F)[opt]

##########  K-M survival curves for MNB in TCGA and CGGA datasets
groups=matrix(0,1,length(status0))
groups[predicted0<=sort(predicted0,F)[opt]]=1
groups[predicted0>sort(predicted0,F)[opt]]=2
# groups[predicted0<=median(predicted0)]=1
# groups[predicted0>median(predicted0)]=2
groups=t(groups)
groups=as.numeric(groups)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


fit<- survfit(Surv(time, status) ~ groups, data = as.data.frame(Survival_Gene))

ggsurv <- ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
                     palette = c("#68BFCF","#F94040"),
                     xlab='Time (days)', xlim=c(0,5000),break.x.by=1000, 
                     ylim=c(0,1), break.y.by=0.5,
                     font.x=c(12,"bold"), font.y=c(12,"bold"), font.tickslab=c(12, "bold"),
                     font.main=c(13,"bold"),
                     risk.table = F, legend=c(0.5,0.95),
                     legend.title=c(""),
                     legend.labs=c("Tim-3 Low (n=403)","Tim-3 High (n=260)"),
                     pval.coord=c(150,0.25))
ggsurv


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


my_comparisons <- list( c("GBM", "Normal"), c("LGG", "Normal") )
ggviolin(Vio_data, x = "group", y = "expression", fill = "group",
         palette = c("#FC4E07","#E7B800","#00AFBB"),
         add = "boxplot", add.params = list(fill = "white"), facet.by = "gene",
         xlab = "(n.GBM=152, n.LGG=511, n.Normal=207)",ylab = "Gene Expression",
         font.x=c(12,"bold"), font.y=c(12,"bold"), font.tickslab=c(11, "bold"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")  # Add significance levels






