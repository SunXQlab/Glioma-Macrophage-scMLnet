###########################################################################
#
# integration 
#
###########################################################################

## library 

library(Seurat)
library(dplyr)

library(patchwork)
library(ggplot2)

library(ggsci)

rm(list = ls())
gc()

## color 

scales::show_col(pal_npg(palette = "nrc", alpha = 0.5)(10))
mycolors_npg <- pal_npg(palette = "nrc", alpha = 0.5)(10)

mycolor_clu <- mycolors_npg
names(mycolor_clu) <- c("1","0","3","5","2","6","7","4","8")

mycolor_ct <- mycolors_npg[c(1:4,6,7,9,10)]
names(mycolor_ct) <- c("Malignant","Macrophages","Fibroblasts","Tcells",
                        "Neutrophil","Endothelial","Oligodendrocytes")
scales::show_col(mycolor_ct)

## function 

StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


## creat seur.obj 

dirs <- list.dirs('./data')[-1]
pten.seur.list <- lapply(dirs,function(dir){
  pten.mtx <- ReadMtx(
    mtx = paste0(dir,"/matrix.mtx.gz"), 
    features = paste0(dir,"/features.tsv.gz"),
    cells = paste0(dir,"/barcodes.tsv.gz"),
  )
  pten.seur <- CreateSeuratObject(counts = pten.mtx,project = gsub('.*/|_.*','',dir))
})
names(pten.seur.list) <- gsub('.*/|_.*','',dirs)

## Standard pre-processing workflow 

pten.list <- list()
for (patient in names(pten.seur.list)) {
  
  pten.patient <- pten.seur.list[[patient]]
  head(pten.patient@meta.data)
  
  pten.patient[["percent.mt"]] <- PercentageFeatureSet(pten.patient, pattern = "^MT-")
  VlnPlot(pten.patient, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  pten.patient <- subset(pten.patient, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
  
  pten.patient <- pten.patient %>%
    SCTransform(vars.to.regress = "percent.mt") %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters() %>% 
    FindClusters(resolution = 0.2)
  # head(pten.patient@meta.data)
  
  pten.list[[patient]] <- pten.patient
  
}

## integration 

features <- SelectIntegrationFeatures(object.list = pten.list, nfeatures = 3000)
pten.list <- PrepSCTIntegration(object.list = pten.list, anchor.features = features)

pten.anchors <- FindIntegrationAnchors(object.list = pten.list, 
                                       normalization.method = "SCT", 
                                       anchor.features = features)
pten.combined <- IntegrateData(anchorset = pten.anchors, normalization.method = "SCT")

DefaultAssay(pten.combined) <- "integrated"
pten.combined <- RunPCA(pten.combined, verbose = FALSE, npcs = 30)
pten.combined <- RunUMAP(pten.combined, reduction = "pca", dims = 1:30)
pten.combined <- FindNeighbors(pten.combined, reduction = "pca", dims = 1:30)
pten.combined <- FindClusters(pten.combined,resolution = 0.1)

p1 <- DimPlot(pten.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(pten.combined, reduction = "umap", group.by = "integrated_snn_res.0.1", label = TRUE,
              repel = TRUE)
p1 + p2

## annotation 

DefaultAssay(pten.combined) <- 'SCT'
Idents(pten.combined) <- pten.combined$integrated_snn_res.0.1

classmk <- c('SOX2','OLIG2', 'ASCL1', # Tumor
             'CD3D','CD3E', # Tcell
             'MAG', 'MOG', # Oligodendrocytes
             'CD14', 'CD68', # Macrophages
             'CSF1R', 'FCGR3A', # Macrophages
             'S100A9','S100A8', # Neutrophil 
             "FBLN1","DCN", # Fibroblasts
             "CLDN5",'VWF' # Endothelial
)
pt.marker <- StackedVlnPlot(pten.combined, c(classmk), pt.size=0, cols=my36colors)
pt.marker
ggsave('./data/pten.combined.cluster.markers.pdf',width = 8,height = 20)

Clu2CT <- function(label){
  
  ct <- switch(label,
               '0' = 'Macrophages',
               '1' = 'Malignant cells',
               '2' = 'Malignant cells',
               '3' = 'Fibroblasts',
               '4' = 'Malignant cells',
               '5' = 'Tcells',
               '6' = 'Neutrophil',
               '7' = 'Endothelial cells',
               '8' = 'Oligodendrocytes'
  )
  
  return(ct)
  
}
pten.combined$celltype <- lapply(as.character(pten.combined$integrated_snn_res.0.1),Clu2CT) %>% unlist()

p1 <- DimPlot(pten.combined, reduction = "umap", group.by = "orig.ident", cols=mycolors_npg[1:4]) +
  labs(title = 'patient')
p2 <- DimPlot(pten.combined, reduction = "umap", group.by = "celltype", label = TRUE,
              repel = TRUE, cols=mycolor_ct, label.box = T) + labs(title = 'cell type')
p1 + p2

Idents(pten.combined) <- pten.combined$celltype
pt.marker <- StackedVlnPlot(pten.combined, c(classmk), pt.size=0, cols=my36colors)
pt.marker
ggsave('./data/pten.combined.celltype.markers.pdf',width = 8,height = 20)

pd <- pten.combined@meta.data
write.table(pd,'./data/pten.combined.final.metadata.txt',sep = '\t',quote = F)

saveRDS(pten.combined,file = './data/pten.final.rds')

###########################################################################
#
# CNV analysis 
#
###########################################################################

## library 

library(rjags)
library(infercnv)

library(Seurat)
library(dplyr)

library(patchwork)
library(ggplot2)

library(AnnoProbe)
library(ggsci)

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

rm(list = ls())
gc()

## color 

scales::show_col(pal_npg(palette = "nrc", alpha = 0.5)(10))
mycolors_npg <- pal_npg(palette = "nrc", alpha = 0.5)(10)

mycolor_ct <- mycolors_npg[c(1:4,6,7,9,10)]
names(mycolor_ct) <- c("Malignant","Macrophages","Fibroblasts","Tcells",
                       "Neutrophil","Endothelial","Oligodendrocytes")
scales::show_col(mycolor_ct)

mycolor_clu <- mycolors_npg
names(mycolor_clu) <- c("1","0","3","5","2","6","7","4","8")

## input 

pten.final <- readRDS(file = './data/pten.final.rds')
DefaultAssay(pten.final) <- 'SCT'
table(pten.final$final_label)
# Ref: Fibroblasts + Tcells
# Que: Malignant cells

malignant.cells  <- row.names(pten.final@meta.data)[which(pten.final$final_label=='Glial cells')]
length(malignant.cells)
malignant.Mat <- as.data.frame(GetAssayData(subset(pten.final, cells=malignant.cells)))

fib.cells <- row.names(pten.final@meta.data)[which(pten.final$final_label=='Fibroblasts')]
fib.cells <- sample(fib.cells,800)
fibMat <- as.data.frame(GetAssayData(subset(pten.final, cells=fib.cells)))

t.cells <- row.names(pten.final@meta.data)[pten.final$final_label=='Tcells']
t.cells <- sample(t.cells,800)
tMat <- as.data.frame(GetAssayData(subset(pten.final, cells=t.cells)))

dat <- cbind(malignant.Mat,fibMat,tMat)
groupinfo <- data.frame(v1=colnames(dat),
                        v2=c(rep('query-malignant',ncol(malignant.Mat)),
                             rep('spike-fib',300),
                             rep('refer-fib',500),
                             rep('spike-t',300),
                             rep('refer-t',500)))

geneInfor <- annoGene(rownames(dat),"SYMBOL",'human')
geneInfor <- geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

dat <- dat[rownames(dat) %in% geneInfor[,1],]
dat <- dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)

write.table(dat,file = './cnv/expFile.txt',sep = '\t',quote = F)
write.table(groupinfo,file = './cnv/groupFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file = './cnv/geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)

## run

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = './cnv/expFile.txt',
                                    gene_order_file = './cnv/geneFile.txt',
                                    annotations_file = './cnv/groupFiles.txt',
                                    delim = "\t", ref_group_names = c("refer-t","refer-fib")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir='./cnv/result/', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

saveRDS(infercnv_obj, './cnv/infercnv_obj.rds')

## figure

infercnv_obj = readRDS("./cnv/result/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data

normal_loc <- infercnv_obj@reference_grouped_cell_indices
test_loc <- infercnv_obj@observation_grouped_cell_indices

anno.df = data.frame(
  CB = c(colnames(expr)[unlist(normal_loc)],colnames(expr)[unlist(test_loc)]),
  class = c(rep("refer-t",length(normal_loc$`refer-t`)),
            rep("refer-fib",length(normal_loc$`refer-fib`)),
            rep("query-malignant",length(test_loc$`query-malignant`)),
            rep("spike-fib",length(test_loc$`spike-fib`)),
            rep("spike-t",length(test_loc$`spike-t`)))
)
head(anno.df)

gn <- rownames(expr)
geneFile <- read.table("./cnv/geneFile.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]

## kmeans clustering

set.seed(20211219)
kmeans.result <- kmeans(t(expr), 4)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB <- rownames(kmeans_df)
kmeans_df <- kmeans_df%>%inner_join(anno.df,by="CB") 
kmeans_df_s <- arrange(kmeans_df,kmeans_class) 
rownames(kmeans_df_s) <- kmeans_df_s$CB
kmeans_df_s$CB <- NULL
kmeans_df_s$kmeans_class <- as.factor(kmeans_df_s$kmeans_class) 
head(kmeans_df_s)
kmeans_df_s <- kmeans_df_s[order(kmeans_df_s$class,kmeans_df_s$kmeans_class),]
kmeans_df_s <- kmeans_df_s[,c(2,1)]

top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v <- RColorBrewer::brewer.pal(8, "Dark2")[1:4]
names(color_v) <- as.character(1:4)
color_c <- RColorBrewer::brewer.pal(8, "Set2")[1:5]
names(color_c) <- unique(anno.df$class)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=color_c,kmeans_class=color_v))

pdf("./cnv/res_kmeans_celltype.pdf",width = 10,height = 6)
ht = Heatmap(matrix = t(expr)[rownames(kmeans_df_s),],
             col = colorRamp2(c(0.9,1,1.1), c("#377EB8","#F0F0F0","#E41A1C")), 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")),
             column_gap = unit(2, "mm"),
             heatmap_legend_param = list(title = "Modified expression",
                                         direction = "vertical",
                                         title_position = "leftcenter-rot",
                                         at=c(0.9,1,1.1),
                                         legend_height = unit(3, "cm")),
             top_annotation = top_anno,left_annotation = left_anno, 
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

## cnv score

expr2 <- expr-1
expr2 <- expr2 ^ 2
CNV_score <- as.data.frame(colMeans(expr2))
colnames(CNV_score) <- "CNV_score"
CNV_score$CB <- rownames(CNV_score)
kmeans_df_s$CB <- rownames(kmeans_df_s)
CNV_score <- CNV_score%>%inner_join(kmeans_df_s,by="CB")

p1 <- ggplot(CNV_score,aes(class,CNV_score))+
  geom_violin(aes(fill=class),color="NA")+
  scale_fill_manual(values = color_c)+
  theme_classic() + ylim(c(0,0.01)) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = 'none'
  )

CNV_score$group <- 'normal cell'
CNV_score$group[CNV_score$class == 'query-malignant'] <- 'malignant cell'
color_g <- color_c[2:3]
names(color_g) <- c('normal cell','malignant cell')

p2 <- ggplot(CNV_score,aes(group,CNV_score))+
  geom_violin(aes(fill=group),color="NA")+
  scale_fill_manual(values = color_g)+
  theme_classic() + ylim(c(0,0.01)) +
  ggsignif::geom_signif(comparisons = list(c("normal cell", "malignant cell")),
                        map_signif_level=T,
                        textsize=6,test=t.test,step_increase=0.2) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 0.5),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p1 + p2 + plot_layout(ncol = 2, widths = c(1, 0.5))
ggsave("./cnv/res_kmeans_celltype_score.pdf",width = 5,height = 3.8)

###########################################################################
#
# CCC Inference 
#
###########################################################################

## library

library(Seurat)
library(dplyr)

library(patchwork)
library(SeuratWrappers)

library(OmnipathR)

library(Matrix)
library(dplyr)

library(ggsci)
library(RColorBrewer)
library(scales)

library(ggplot2)
library(plotrix)
library(igraph)
library(ggraph)

rm(list = ls())
gc()

## color 

scales::show_col(pal_npg(palette = "nrc", alpha = 0.5)(10))
mycolors_npg <- pal_npg(palette = "nrc", alpha = 0.5)(10)

mycolor_ct <- mycolors_npg[c(1:4,6,7,9,10)]
names(mycolor_ct) <- c("Malignant","Macrophages","Fibroblasts","Tcells",
                       "Neutrophil","Endothelial","Oligodendrocytes")
scales::show_col(mycolor_ct)

mycolor_key <- brewer.pal(4, "Set2")
names(mycolor_key) <- c('Ligand','Receptor','Target','TF')
scales::show_col(mycolor_key)

## function 

fisher_test <- function(subset1,subset2,backgrond)
{
  a=length(intersect(subset1,subset2))
  b=length(subset1)-a
  c=length(subset2)-a
  d=length(backgrond)-a-b-c
  matrix=matrix(c(a,c,b,d),nrow=2)
  fisher.test(matrix,alternative="greater")$p.value
}

getTFTG <- function(TFTG.DB, target.degs, target.genes)
{
  
  if (!is.data.frame(TFTG.DB))
    stop("TFTG.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(TFTG.DB))
    stop("TFTG.DB must contain a column named 'source'")
  if (!"target" %in% colnames(TFTG.DB))
    stop("TFTG.DB must contain a column named 'target'")
  
  # get TF list
  TF.list <- TFTG.DB %>% select(source) %>% unlist() %>% unique()
  
  # get Target list
  TG.list <- lapply(TF.list, function(x){
    TFTG.DB %>% filter(source == x)  %>% select(target) %>% unlist() %>% unique() %>% intersect(.,target.genes)
  })
  names(TG.list) <- TF.list
  
  # get target differently expressed genes
  DEGs <- target.degs
  
  # perform fisher test
  TFs <- lapply(TG.list, function(x){fisher_test(subset1 = x, subset2 = DEGs, backgrond = target.genes)})
  TFs <- unlist(TFs)
  TFs <- names(TFs)[TFs <= 0.05]
  TFs <- TFs[TFs %in% target.genes]
  
  # get activated LR pairs
  TFTGList <- TG.list[TFs]
  TFTGList <- lapply(TFTGList, function(x){intersect(x, DEGs)})
  TFTGList <- paste(rep(TFs, times = lengths(TFTGList)), unlist(TFTGList), sep = "_")
  
  # check result
  if(length(TFTGList)==0)
    stop("Error: No significant TFTG pairs")
  
  # get result
  TFTGTable <- TFTGList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(TFTGTable) <- c("source","target")
  
  cat(paste0("get ",length(TFTGList)," activated TFTG pairs\n"))
  return(TFTGTable)
  
}

getRecTF <- function(RecTF.DB, Rec.list, TF.list)
{
  
  if (!is.data.frame(RecTF.DB))
    stop("RecTF.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'source'")
  if (!"target" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'target'")
  
  # make sure Rec.list in RecTF.DB
  Rec.list <- Rec.list[Rec.list %in% RecTF.DB$source]
  Rec.list <- as.vector(Rec.list)
  
  # make sure TF.list in RecTF.DB
  TF.list <- TF.list[TF.list %in% RecTF.DB$target]
  TF.list <- as.vector(TF.list)
  
  # get TF activated by Receptors
  TFofRec <- lapply(Rec.list, function(x){
    RecTF.DB %>% filter(source == x)  %>% select(target) %>% unlist() %>% unique()
  })
  names(TFofRec) <- Rec.list
  
  # get all TF
  TFofALL <- RecTF.DB %>% select(target) %>% unlist() %>% unique()
  
  # perform fisher test
  Recs <- lapply(TFofRec, function(x){
    fisher_test(subset1 = x, subset2 = TF.list, backgrond = TFofALL)
  })
  Recs <- unlist(Recs)
  Recs <- names(Recs)[Recs <= 0.05]
  
  # get activated RecTF pairs
  RecTFList <- TFofRec[Recs]
  RecTFList <- lapply(RecTFList, function(x){intersect(x, TF.list)})
  RecTFList <- paste(rep(Recs, times = lengths(RecTFList)), unlist(RecTFList), sep = "_")
  
  # check result
  if(length(RecTFList)==0)
    stop("Error: No significant RecTF pairs")
  
  # get result
  RecTFTable <- RecTFList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(RecTFTable) <- c("source","target")
  
  cat(paste0("get ",length(RecTFList)," activated RecTF pairs\n"))
  return(RecTFTable)
  
}

runMLnet <- function(sender,reciver,meanExprMat,rank.cutoff=0.2)
{
  
  ## get LigRec #################
  
  LigRec <- LigRecDB
  LigRec$LRscore <- apply(LigRec,1,function(x){
    
    meanExprMat[x[1],sender]*meanExprMat[x[2],reciver]
    
  })
  LigRec$rank <- rank(-LigRec$LRscore) 
  LigRec <- LigRec[LigRec$rank<quantile(LigRec$rank,rank.cutoff),]
  LigRec <- LigRec[LigRec$target %in% DEGs_list[[reciver]],]
  
  ## get TFTG #################
  
  target_gene <- rownames(meanExprMat)
  target_deg <- DEGs_list[[reciver]]
  tryCatch({
    TFTG <- getTFTG(TFTGDB, target_deg, target_gene)
  },
  error=function(e){cat(conditionMessage(e),"\n")})
  tag1 = exists("TFTG")
  if(!tag1) return(NULL)
  
  ## get RecTF ################
  
  Rec.list <- LigRec$target %>% as.character() %>% unique()
  TF.list <- TFTG$source %>% as.character() %>% unique()
  tryCatch({
    RecTF <- getRecTF(RecTFDB, Rec.list, TF.list)
  },
  error=function(e){cat(conditionMessage(e),"\n")})
  tag2 = exists("RecTF")
  if(!tag2) return(NULL)
  
  ## connect ##################
  
  Receptors <- intersect(LigRec$target, RecTF$source)
  TFs <- intersect(RecTF$target, TFTG$source)
  
  LigRec <- LigRec[LigRec$target %in% Receptors,] %>% dplyr::select(source, target)
  RecTF <- RecTF[RecTF$source %in% Receptors & RecTF$target %in% TFs,] %>% dplyr::select(source, target)
  TFTG <- TFTG[TFTG$source %in% TFs,] %>% dplyr::select(source, target)
  
  ## multi-layer ############
  
  LigRec <- LigRec[!duplicated(LigRec),]
  RecTF <- RecTF[!duplicated(RecTF),]
  TFTG <- TFTG[!duplicated(TFTG),]
  
  MLnet <- list(LigRec = LigRec,
                RecTF = RecTF,
                TFTG = TFTG)
  
  return(MLnet)
  
}

DrawLeg <- function(alltype,allcolos,xstart,ystart,cirr,jiange)
{
  ThisX <- xstart
  ThisY <- ystart
  for(i in seq(1,length(alltype)))
  {
    ThisType <- alltype[i]
    draw.circle(ThisX,ThisY,cirr,col = allcolos[ThisType])
    text(ThisX+cirr+0.1,ThisY,ThisType,adj = 0)
    
    ThisY <- ThisY - (2*cirr + jiange)
  }
}

DrawCellComm <- function(CellTab,colodb)
{
  
  aaa <- CellTab
  g <- graph.data.frame(aaa,directed=TRUE)
  
  alltype <- c(aaa$cell_from,aaa$cell_to)
  alltype <- unique(alltype)
  if(is.null(names(colodb))){
    allcolos <- colodb
    names(allcolos) <- alltype
  }else{
    allcolos <- colodb
  }
  
  edge.start <- ends(g,es=E(g),names=FALSE)
  
  layout <- in_circle()  #igraph packages
  coords <- layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  
  loop.angle <- ifelse(coords_scale[V(g),1]>0,-atan(coords_scale[V(g),2]/coords_scale[V(g),1]),pi-atan(coords_scale[V(g),2]/coords_scale[V(g),1]))
  
  vertex.label.color <- 'black'
  V(g)$size <- 20
  V(g)$color <- allcolos[V(g)]
  V(g)$label.color <- vertex.label.color
  V(g)$label.cex <- 0
  
  label <- FALSE
  if(label){
    E(g)$label<-E(g)$n
  }
  
  edge.max.width = 10
  if(max(E(g)$n)==min(E(g)$n)){
    E(g)$width <- 2
  }else{
    E(g)$width <- 1+edge.max.width/(max(E(g)$n)-min(E(g)$n))*(E(g)$n-min(E(g)$n))
  }
  
  edge.label.color <- "black"
  E(g)$arrow.width<- 1
  E(g)$arrow.size <- 0.1
  E(g)$label.color <- edge.label.color
  E(g)$label.cex <- 1
  E(g)$color <- V(g)$color[edge.start[,1]]
  
  if(sum(edge.start[,2]==edge.start[,1])!=0)
  {
    E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  
  #draw 
  shape = 'circle'
  margin = 0.2
  edgecurved = 0.2
  vertexlabelcex = 1
  plot(g,edge.curved=edgecurved,vertex.label = "",vertex.shape=shape,layout=coords_scale,margin=margin,vertex.shape="fcircle",vertex.label.cex = vertexlabelcex,axes = FALSE)
  
  #draw legend
  xstart <- 1.5
  ystart <- 0.5
  jiange <- 0.02
  cirr <- 0.04
  DrawLeg(alltype,allcolos,xstart,ystart,cirr,jiange)
  
}

DrawMLnetwork <- function(df_edges,df_nodes,colodb)
{
  
  subnet <- graph_from_data_frame(d = df_edges, vertices = df_nodes)
  df_nodes <- df_nodes[match(names(V(subnet)),df_nodes$node),]
  root_index <- grep('Ligand',df_nodes$key)
  set.seed(4)
  coords <- layout_(subnet,layout = as_tree(root = root_index)) %>% as.data.frame()
  coords <- cbind(coords,df_nodes$key)
  colnames(coords) <- c('dim_x','dim_y','type')
  
  dist_TG <- 1;len_TG <- table(coords$type)[['Target']]
  dist_TF <- 1.5;len_TF <- table(coords$type)[['TF']]
  dist_Rec <- 2.5;len_Rec <- table(coords$type)[['Receptor']]
  dist_Lig <- 2.5;len_Lig <- table(coords$type)[['Ligand']]
  if(len_Lig==1){
    coords$dim_x[coords$type == 'Ligand'] = 0
  }else{
    dim_x_1 = seq(to = -dist_Lig/2,by = dist_Lig,length.out = ceiling(len_Lig/2))
    dim_x_2 = seq(from = dist_Lig/2,by = dist_Lig,length.out = len_Lig-ceiling(len_Lig/2))
    coords$dim_x[coords$type == 'Ligand'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_Rec==1){
    coords$dim_x[coords$type == 'Receptor'] = 0
  }else{
    dim_x_1 = seq(to = -dist_Rec/2,by = dist_Rec,length.out = ceiling(len_Rec/2))
    dim_x_2 = seq(from = dist_Rec/2,by = dist_Rec,length.out = len_Rec-ceiling(len_Rec/2))
    coords$dim_x[coords$type == 'Receptor'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_TG<len_TF & len_TF>10){
    
    dist_TG <- 1
    dist_TF <- 0.5
    
  }else if(len_TG>len_TF & len_TG>10){
    
    dist_TG <- 0.5
    dist_TF <- 1
    
  }
  if(len_TF==1){
    coords$dim_x[coords$type == 'TF'] = 0
  }else{
    dim_x_1 = seq(to = -dist_TF/2,by = dist_TF,length.out = ceiling(len_TF/2))
    dim_x_2 = seq(from = dist_TF/2,by = dist_TF,length.out = len_TF-ceiling(len_TF/2))
    coords$dim_x[coords$type == 'TF'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_TG==1){
    coords$dim_x[coords$type == 'Target'] = 0
  }else{
    dim_x_1 = seq(to = -dist_TG/2,by = dist_TG,length.out = ceiling(len_TG/2))
    dim_x_2 = seq(from = dist_TG/2,by = dist_TG,length.out = len_TG-ceiling(len_TG/2))
    coords$dim_x[coords$type == 'Target'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  coords$dim_y <- lapply(coords$type,switch,'Ligand'=0.9,'Receptor'=0.6,'TF'=0.3,'Target'=0) %>% unlist()
  # plot(subnet, layout = as.matrix(coords[,1:2]), vertex.color=df_nodes$color)
  
  layout <- create_layout(subnet, layout = 'tree')
  layout[,1:2] <- coords[,1:2]
  # head(layout)
  
  temp <- function(key){
    
    max(coords$dim_x[coords$type==key])+0.5
    
  }
  df_anno <- data.frame(x = lapply(unique(layout$key), temp) %>% unlist(),
                        y = c(0.92,0.62,0.32,0.02),
                        lab = unique(layout$key),
                        key = unique(layout$key))
  
  pt <- ggraph(layout) +
    geom_edge_link(aes(),color="grey",
                   arrow = arrow(length = unit(1.5, 'mm')), 
                   start_cap = circle(3, 'mm'),
                   end_cap = circle(3, 'mm')) + 
    geom_node_point(aes(fill = key,color = key),shape=21,size = 8) +
    geom_node_text(aes(label=name),size=2) + # ,fontface='bold',family='Arial'
    xlim(c(min(layout$x)-2,max(layout$x)+2)) +
    geom_text(data = df_anno, aes(x,y,label=lab,color=key),
              vjust = 1, hjust = 0, size = 4,fontface='bold') + # ,family='ARL'
    scale_fill_manual(values = colodb) + scale_color_manual(values = colodb) +
    guides(fill='none',color='none') +
    theme_graph()
  print(pt)
  
  return(pt)
  
}

## Database 

ppi_form_kinaseextra <- import_kinaseextra_interactions() 
ppi_form_omnipath <- import_omnipath_interactions() 
ppi_form_pathwayextra <- import_pathwayextra_interactions() 
ppi_from_Omnipath <- rbind(ppi_form_kinaseextra,ppi_form_omnipath,ppi_form_pathwayextra) 
ppiDB <- ppi_from_Omnipath %>% dplyr::filter(n_resources>0) %>% 
  dplyr::select(source_genesymbol, target_genesymbol, n_resources) %>%
  dplyr::rename(source=source_genesymbol, target=target_genesymbol, weight = n_resources) %>%
  dplyr::filter(source %in% rownames(meanExprMat), target %in% rownames(meanExprMat))

LigRecDB <- import_ligrecextra_interactions()
LigRecDB <- LigRecDB %>% 
  dplyr::select(source_genesymbol, target_genesymbol, n_resources) %>%
  dplyr::rename(source=source_genesymbol, target=target_genesymbol, weight = n_resources) %>%
  dplyr::filter(source %in% unique(c(ppiDB$source,ppiDB$target)), 
                target %in% unique(c(ppiDB$source,ppiDB$target)))

TFTG_from_dorothea <- import_dorothea_interactions()
TFTG_from_tftarget <- import_tf_target_interactions() 
TFTGDB <- rbind(TFTG_from_dorothea[,-grep("dorothea_level",colnames(TFTG_from_dorothea))],
                TFTG_from_tftarget) 
TFTGDB <- TFTGDB %>% 
  dplyr::select(source_genesymbol, target_genesymbol, n_resources) %>%
  dplyr::rename(source=source_genesymbol, target=target_genesymbol, weight = n_resources) %>%
  dplyr::filter(source %in% unique(c(ppiDB$source,ppiDB$target)), 
                target %in% unique(c(ppiDB$source,ppiDB$target))) 

Database <- list(LigRecDB = LigRecDB,
                 ppiDB = ppiDB,
                 TFTGDB = TFTGDB)
saveRDS(Database, "./cci/database.rds")

V = unique(c(ppiDB$source,ppiDB$target))
adjM <- matrix(0,nrow = length(V),ncol = length(V), dimnames = list(V,V))
for(i in 1:nrow(ppiDB)){
  
  adjM[ppiDB$source[i],ppiDB$target[i]] <- 1/ppiDB$weight[i]
  adjM[ppiDB$target[i],ppiDB$source[i]] <- 1/ppiDB$weight[i]
  
}
adjM <- as(adjM, "dgCMatrix")

G <- igraph::graph_from_adjacency_matrix(adjM, mode='undirected', weighted=TRUE)
G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
isolates <- which(degree(G, mode = c("all")) == 0) - 1
G <- delete.vertices(G, names(isolates))

Receptors <- LigRecDB$target %>% as.character() %>% unique()
TFs <- TFTGDB$source %>% as.character() %>% unique()
distM <- lapply(Receptors, function(rec){
  
  igraph::distances(G, v = rec, to = TFs)
  
}) %>% do.call("rbind",.) %>% as.data.frame()
rownames(distM) <- Receptors
colnames(distM) <- TFs

score.cutoff = 0.2
RecTFDB <- distM
RecTFDB <- cbind(rownames(RecTFDB),RecTFDB)
RecTFDB <- reshape2::melt(RecTFDB, id = "rownames(RecTFDB)")
colnames(RecTFDB) <- c('source','target','score')
RecTFDB <- RecTFDB %>%
  dplyr::filter(score != Inf, score != 0) %>% 
  dplyr::filter(score <= quantile(score, score.cutoff))
saveRDS(RecTFDB, "./cci/RecTFDB.rds")

## input 

seur <- readRDS("./data/pten.final.rds")
metadata <- seur@meta.data
Idents(seur) <- seur$celltype

if(T){
  
  seur_patient <- subset(seur,subset = orig.ident == 'M47')
  metadata_patient <- metadata[metadata$orig.ident == 'M47',]
  
  metadata_patient <- data.frame(barcode = rownames(metadata_patient),
                                 celltype = metadata_patient$celltype)
  cts <- metadata_patient$celltype %>% unique()
  
  meanExprMat <- lapply(cts, function(ct){
    
    exprMat <- seur[['SCT']]@data
    meanExpr <- rowMeans(exprMat[, seur$celltype == ct])
    
  }) %>% do.call("cbind",.) %>% as.data.frame()
  colnames(meanExprMat) <- cts
  
  all.markers <- FindAllMarkers(seur_patient, only.pos = TRUE) 
  DEGs_list <- split(all.markers$gene,all.markers$cluster)
  
  save(cts,meanExprMat,all.markers,DEGs_list,file = paste0('./cci/PTEN_loss/input.rda'))
  
}
if(T){
  
  seur_patient <- subset(seur,subset = orig.ident %in% c('M48','M58'))
  metadata_patient <- metadata[metadata$orig.ident %in% c('M48','M58'),]
  
  metadata_patient <- data.frame(barcode = rownames(metadata_patient),
                                 celltype = metadata_patient$celltype)
  cts <- metadata_patient$celltype %>% unique()
  
  meanExprMat <- lapply(cts, function(ct){
    
    exprMat <- seur[['SCT']]@data
    meanExpr <- rowMeans(exprMat[, seur$celltype == ct])
    
  }) %>% do.call("cbind",.) %>% as.data.frame()
  colnames(meanExprMat) <- cts
  
  all.markers <- FindAllMarkers(seur_patient, only.pos = TRUE) 
  DEGs_list <- split(all.markers$gene,all.markers$cluster)
  
  save(cts,meanExprMat,all.markers,DEGs_list,file = paste0('./cci/WT/input.rda'))
  
}

## run 

groups <- c('PTEN-loss','WT')
for (group in groups) {
  
  cat(paste0(group,'\n'))
  
  ## load
  
  load(paste0("./cci/",group,"/input.rda"))
  
  Database <- readRDS("./cci/database.rds")
  LigRecDB <- Database$LigRecDB
  TFTGDB <- Database$TFTGDB
  RecTFDB <- readRDS("./cci/RecTFDB.rds")
  
  ## group
  
  df_cellpair <- paste(rep(cts,each=length(cts)), rep(cts,times=length(cts)), sep = "_")
  df_cellpair <- lapply(df_cellpair, function(x){strsplit(x,"_")[[1]]}) %>% do.call("rbind",.)
  df_cellpair <- df_cellpair[apply(df_cellpair, 1, function(x){x[1]!=x[2]}),]
  df_cellpair <- data.frame(df_cellpair)
  colnames(df_cellpair) <- c("sender","receiver")
  
  ## run
  MLnet <- apply(df_cellpair, 1, function(cp){
    
    cat(paste0("sender: ",cp[1],'\n'))
    cat(paste0("receiver: ",cp[2],'\n'))
    res <- runMLnet(cp[1],cp[2],meanExprMat,rank.cutoff = 0.2)
    cat("###########################\n")
    res
    
  })
  names(MLnet) <- paste(df_cellpair$sender, df_cellpair$receiver, sep = "_")
  str(MLnet)
  saveRDS(MLnet, paste0("./cci/",group,"/mlnet/MLnet.rds"))
  
}

## cell_comm

colodb <- mycolor_ct
show_col(colodb)

MLnet <- readRDS("./cci/PTEN-loss/mlnet/MLnet.rds")
df_cellpair <- lapply(1:length(MLnet),function(i){
  
  net <- MLnet[[i]]
  dat <- data.frame(cell_from = lapply(names(MLnet)[i], function(x){strsplit(x,"_")[[1]][1]}) %>% unlist(),
                    cell_to = lapply(names(MLnet)[i], function(x){strsplit(x,"_")[[1]][2]}) %>% unlist(),
                    n_TGs = ifelse(is.null(net),0,length(unique(net$TFTG$target))))
  dat
  
}) %>% do.call('rbind',.) %>% as.data.frame()

LRtab <- df_cellpair
colnames(LRtab) <- c('cell_from','cell_to','n')

pdf(paste0("./cci/M47/figure/cell_comm.pdf"),height = 5,width = 6)
DrawCellComm(LRtab,colodb)
dev.off()

## Multilayer network

mlnet <- MLnet$Malignant_Macrophages
xlsx::write.xlsx(mlnet$LigRec,file = './cci/M47/mlnet/MLnet_Malignant_Macrophages.xlsx',sheetName = 'LigRec',row.names = F)
xlsx::write.xlsx(mlnet$RecTF,file = './cci/M47/mlnet/MLnet_Malignant_Macrophages.xlsx',sheetName = 'RecTF',row.names = F,append = T)
xlsx::write.xlsx(mlnet$TFTG,file = './cci/M47/mlnet/MLnet_Malignant_Macrophages.xlsx',sheetName = 'TFTG',row.names = F,append = T)

Database <- readRDS("./cci/database.rds")
LigRecDB <- Database$LigRecDB
lig = "LGALS9"
rec = LigRecDB$target[LigRecDB$source == lig]

mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$source %in% lig & mlnet$LigRec$target %in% rec,]
mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$source %in% mlnet$LigRec$target,]
mlnet$TFTG <- mlnet$TFTG[mlnet$TFTG$source %in% mlnet$RecTF$target,]
xlsx::write.xlsx(mlnet$LigRec,file = './cci/M47/mlnet/MLnet_LGALS9_HAVCR2.xlsx',sheetName = 'LigRec',row.names = F)
xlsx::write.xlsx(mlnet$RecTF,file = './cci/M47/mlnet/MLnet_LGALS9_HAVCR2.xlsx',sheetName = 'RecTF',row.names = F,append = T)
xlsx::write.xlsx(mlnet$TFTG,file = './cci/M47/mlnet/MLnet_LGALS9_HAVCR2.xlsx',sheetName = 'TFTG',row.names = F,append = T)

m2m <- read.csv('./cci/M2 marker genes.csv',header = F)
m2m <- unlist(m2m) %>% .[.!='']
human <- biomaRt::useMart('ensembl',dataset = 'hsapiens_gene_ensembl')
mouse <- biomaRt::useMart('ensembl',dataset = 'mmusculus_gene_ensembl')
m2m <- biomaRt::getLDS(attributes = 'mgi_symbol',filters = 'mgi_symbol',
                       values = m2m, mart = mouse,
                       attributesL = 'hgnc_symbol',martL = human,uniqueRows = T)
m2m <- m2m$HGNC.symbol

mlnet <- MLnet$Malignant_Macrophages
mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$source %in% lig & mlnet$LigRec$target %in% rec,]
mlnet$TFTG <- mlnet$TFTG[mlnet$TFTG$target %in% m2m,]
mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$source %in% mlnet$LigRec$target &
                             mlnet$RecTF$target %in% mlnet$TFTG$source,]
mlnet$TFTG <- mlnet$TFTG[mlnet$TFTG$source %in% mlnet$RecTF$target,]

df_edges <- do.call("rbind",mlnet)
df_nodes <- data.frame(node = c(unique(mlnet$LigRec$source),
                                unique(mlnet$LigRec$target),
                                unique(mlnet$TFTG$source),
                                unique(mlnet$TFTG$target)))
df_nodes$key <- c(rep('Ligand',length(unique(mlnet$LigRec$source))),
                  rep('Receptor',length(unique(mlnet$LigRec$source))),
                  rep('TF',length(unique(mlnet$TFTG$source))),
                  rep('Target',length(unique(mlnet$TFTG$target))))
mycolor_node <- mycolor_key[df_nodes$key]

pt_network <- DrawMLnetwork(df_edges,df_nodes,mycolor_key)
ggsave(plot = ,'./cci/WT/figure/multilayer_network.pdf')

## check expr

pten.final <- readRDS("./data/pten.final.rds")
DefaultAssay(pten.final) <- 'SCT'
pten.final@meta.data <- metadata

pten.final$stage <- 'WT'
pten.final$stage[pten.final$orig.ident == 'M47'] <- 'PTEN_LOSS'

StackedVlnPlot(pten.final, c(Lig,Recs), pt.size=0, cols=my36colors,group.by = 'stage')
ggsave('./expr/pten.final.vln.stage.pdf',width = 3,height = 13)

## check celltype propotion

DimPlot(pten.final, reduction = "umap", split.by = "stage", group.by = "celltype",
        label = T, repel = T, label.box = T, cols = mycolor_ct2) + NoLegend()
ggsave('./expr/pten.final.celltype.pdf',width = 10,height = 5)

df_cellnum <- as.data.frame(as.matrix(table(pten.final$stage,pten.final$celltype)))
df_cellnum <- df_cellnum %>% dplyr::filter(Var1 == 'PTEN_LOSS') %>%
  mutate(Group = factor(Var2),
         cumulative = cumsum(Freq),
         midpoint = cumulative - Freq / 2,
         label = paste0(Group, " ", round(Freq / sum(Freq) * 100, 1), "%"))
df_cellnum <- df_cellnum[order(df_cellnum$Freq),]

p1 <- ggplot(df_cellnum, aes(x = 1, weight = Freq, fill = Group)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") + ## 以y轴建立极坐标
  scale_fill_manual(values = mycolor_ct2) + 
  geom_text(aes(x = 1.3, y = midpoint, label = label)) + ## 加上百分比标签的位置和数值
  ggmap::theme_nothing()   

pdf('./data/pten.loss.celltype.proportion.pdf',width = 4,height = 4)
p1
dev.off()

df_cellnum <- as.data.frame(as.matrix(table(pten.final$stage,pten.final$celltype)))
df_cellnum <- df_cellnum %>% dplyr::filter(Var1 == 'WT') %>%
  mutate(Group = factor(Var2),
         cumulative = cumsum(Freq),
         midpoint = cumulative - Freq / 2,
         label = paste0(Group, " ", round(Freq / sum(Freq) * 100, 1), "%"))
df_cellnum <- df_cellnum[order(df_cellnum$Freq),]

p2 <- ggplot(df_cellnum, aes(x = 1, weight = Freq, fill = Group)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y") + ## 以y轴建立极坐标
  scale_fill_manual(values = mycolor_ct2) + 
  geom_text(aes(x = 1.3, y = midpoint, label = label)) + ## 加上百分比标签的位置和数值
  ggmap::theme_nothing()   

pdf('./data/pten.wt.celltype.proportion.pdf',width = 4,height = 4)
p2
dev.off()

###########################################################################
#
# Enrichment analysis
#
###########################################################################

## library 

rm(list = ls())
gc()

library(Seurat) 
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

## color

mycols <- c('white','#F8766D')
scales::show_col(mycols)

## genelist

MLnet <- readRDS("./cci/M47/mlnet/MLnet.rds")
MLnet <- MLnet$Malignant_Macrophages
tfs <- MLnet$RecTF$target[MLnet$RecTF$source %in% "HAVCR2"]
tgs <- MLnet$TFTG$target[MLnet$TFTG$source %in% tfs]
geneList <- tgs

g2s <- toTable(org.Hs.egSYMBOL)
geneList <- g2s$gene_id[match(geneList,g2s$symbol)]
geneList <- na.omit(geneList)

## GO

enrich_gobp <- enrichGO(geneList, 'org.Hs.eg.db', ont="BP", keyType = 'ENTREZID', 
                        minGSSize = 1, pvalueCutoff = 0.99)
str(enrich_gobp)

res_enrich_gobp <- enrich_gobp@result
res_enrich_gobp <- res_enrich_gobp[res_enrich_gobp$p.adjust <= 0.05,]

keywords <- 'macrophage activation|macrophage differentiation|macrophage migration|vessel|vasculature'
hits <- res_enrich_gobp$Description[grep(keywords,res_enrich_gobp$Description,ignore.case = T)]
res_enrich_gobp_hits <- res_enrich_gobp[res_enrich_gobp$Description %in% hits,]
res_enrich_gobp_hits$GeneRatio <- res_enrich_gobp_hits$Count/169

pdf("./cci/figure/enrich_gobp.pdf", width = 7, height = 4.5)
p1 <- ggplot(res_enrich_gobp_hits, aes(x = GeneRatio, y = reorder(Description ,GeneRatio), 
                                       size = Count, fill = p.adjust)) + 
  geom_point(shape = 21) + ylab("GO Term") + theme_bw() +
  scale_fill_continuous(low = mycols[2], high = mycols[1]) + 
  theme(axis.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.box = "horizontal") 
p1
dev.off()

###########################################################################
#
# GSEA analysis
#
###########################################################################

## library

rm(list = ls())
gc()

library(Seurat)
library(clusterProfiler)
library(enrichplot)

library(GSEABase) 

library(ggsci)

## color

mycols <- pal_lancet("lanonc", alpha = 0.7)(9)

## genelist

seur <- readRDS("./data/pten.final.rds")
DefaultAssay(seur) <- 'SCT'

seur$stage <- 'WT'
seur$stage[seur$orig.ident == 'M47'] <- 'PTEN_LOSS'

seur_macro <- subset(seur,subset = final_label == 'Macrophages')
Idents(seur_macro) <- seur_macro$stage

degs.macro <- FindMarkers(seur_macro, ident.1 = "PTEN_LOSS", logfc.threshold = 0.01, only.pos = FALSE)
saveRDS(degs.macro,'./cci/pten.macro.markers.lossVSwt.rds')

geneList <-  degs.macro
geneList <- geneList$avg_log2FC
names(geneList) = rownames(degs.macro)
geneList <- geneList[degs.macro$pct.1>=0.1]
geneList <- geneList[order(geneList, decreasing = T)]

## geneset

gmtfile ='./msigDB/c5.go.bp.v7.4.symbols.gmt'
geneset <- read.gmt(gmtfile)
geneset$term <- geneset$term %>% tolower()
length(unique(geneset$term))

gsea_gobp <- GSEA(geneList, TERM2GENE=geneset, minGSSize = 1, pvalueCutoff = 0.99, verbose=FALSE, seed = 10)
str(gsea_gobp) 

res_gsea_gobp <- gsea_gobp@result
res_gsea_gobp <- res_gsea_gobp[res_gsea_gobp$p.adjust <= 0.05 & res_gsea_gobp$setSize > 10,]

hits <- res_gsea_gobp$ID[grep("macrophage",res_gsea_gobp$ID,ignore.case = T)]
pdf("./cci/figure/gsea_lossVSwt_macrophages.pdf", width = 15, height = 10)
gseaplot2(gsea_gobp, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
dev.off()

hits <- res_gsea_gobp$ID[grep("vessel|vascul",res_gsea_gobp$ID,ignore.case = T)]
pdf("./cci/figure/gsea_lossVSwt_vessel_development.pdf", width = 15, height = 10)
gseaplot2(gsea_gobp, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
dev.off()
