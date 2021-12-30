library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratWrappers)

library(ggsci)
library(scales)
mycols <- pal_lancet("lanonc", alpha = 0.7)(9)
show_col(pal_lancet("lanonc", alpha = 0.7)(9))

mycols_scale <- pal_igv("default", alpha = 0.7)(48)
show_col(pal_igv("default", alpha = 0.7)(48))

###############################
## scRNA-seq data processing ##
###############################
## download ##########

dat <- data.table::fread("./2019_cell/downloaded_data/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz")
dim(dat)
dat[1:4,1:4]

dat <- as.data.frame(dat)
rownames(dat) <- dat$GENE
dat <- dat[,-1]

## run seurat ########

wrap_seurat <- function(dat, label){
  
  # creat
  seur.obj <- CreateSeuratObject(counts = dat, project = label)
  
  # quality control (skip)
  # seur.obj[["percent.mt"]] <- PercentageFeatureSet(seur.obj, pattern = "^MT-")
  # seur.obj <- subset(seur.obj, subset = nFeature_RNA > 200 & percent.mt < 10) 
  cat(paste0("cells: ",dim(seur.obj)[2],'\n'))
  cat(paste0("genes: ",dim(seur.obj)[1],'\n'))
  
  # normalize + find HVGs + scale
  seur.obj <- NormalizeData(seur.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seur.obj <- FindVariableFeatures(seur.obj, selection.method = "vst", nfeatures = 1000)
  
  # scale
  
  all.genes <- rownames(seur.obj)
  seur.obj <- ScaleData(seur.obj, features = all.genes)
  
  # Dimensionality Reduction
  
  hvgs <- VariableFeatures(seur.obj)
  seur.obj <- RunPCA(seur.obj, features = hvgs)
  
  # Visulization
  
  seur.obj <- FindNeighbors(seur.obj)
  seur.obj <- RunTSNE(seur.obj)
  
}
seur <- wrap_seurat(dat, "GBM")

## cluster ########

seur <- FindClusters(seur, resolution = 0.1)
DimPlot(seur, reduction = "tsne", pt.size = 1)

## imputation ########

seur <- RunALRA(seur)
DefaultAssay(seur) <- "alra"

## markers ########

## macrophages: CD14, AIF1, FCER1G, FCGR3A, TYROBP, CSF1R
## T cells: CD2, CD3D, CD3E, CD3G
## oligodendrocytes: MBP, TF, PLP1, MAG, MOG, CLDN11
## Malignant: PROM1, MSI1, NES, SOX2, FUT4, PARP1, IFAP

cluster_markers <- c('CD14', 'AIF1', 'FCER1G', 'FCGR3A', 'TYROBP', 'CSF1R', # macrophages
                     'CD2', 'CD3D', 'CD3E', 'CD3G', # T cells
                     'MBP', 'TF', 'PLP1', 'MAG', 'MOG', 'CLDN11', # oligodendrocytes
                     'PROM1', 'MSI1', 'NES', 'SOX2', 'PARP1', 'MBTPS2' # Malignant
                     )
cluster_markers[!cluster_markers %in% rownames(seur[['alra']]@data)]
# pdf(file="./result/marker_heatmap.pdf", width = 10, height = 4)
pp = DotPlot(seur, features = cluster_markers, cols = c('white','#F8766D'), dot.scale =5) + RotatedAxis()
pp = pp + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
    ) + 
  labs(x='',y='') +
  guides(color = guide_colorbar(title = 'Scale expression'),
         size = guide_legend(title = 'Percent expressed')) +
  theme(axis.line = element_line(size = 0.6))
print(pp)
# dev.off()

## annotation ########

seur$orig.annotation <- Idents(seur)
seur$orig.annotation <- factor(seur$orig.annotation, 
                               labels = c('Malignant','Malignant','macrophages','Malignant',
                                          'Malignant','oligodendrocytes','Malignant','Tcell'))
table(seur$orig.annotation, seur$seurat_clusters)

Idents(seur) <- seur$orig.annotation
seur$patient <- colnames(seur) %>% lapply(.,function(x){strsplit(x,"-")[[1]][1]}) %>% unlist()

## plot ########

mycols_seur <- mycols[c(1,2,3,5)]
names(mycols_seur) <- c('Tcell','macrophages','oligodendrocytes','Malignant')
p1 <- DimPlot(seur, reduction = "tsne", pt.size = 1, cols = mycols_seur, 
              label = T, repel = T, label.size = 7, label.box = T, label.color = "white") + 
  NoLegend() + labs(title = "celltype") + theme(plot.title = element_text(hjust = 0.5))
mycols_pat <- mycols_scale[1:length(unique(seur$patient))]
names(mycols_pat) <- unique(seur$patient)
p2 <- DimPlot(seur, reduction = "tsne", pt.size = 1, group.by = "patient", cols = mycols_pat)

# pdf("./result/dimplot_annotation.pdf", width = 16)
p1+p2
# dev.off()

## DEGs ########

All.markers.alra <- FindAllMarkers(seur, only.pos = TRUE) 
write.csv(All.markers.alra, "./result/all.markers.alra.csv")

## save ########
saveRDS(seur, "./result/seur.rds")

#############################
## cell-cell communication ##
#############################

rm(list = ls())

## load input ###############

seur <- readRDS("./result/seur.rds")

metadata <- data.frame(barcode = rownames(seur@meta.data),
                       celltype = seur@active.ident)

cts <- metadata$celltype %>% unique()
meanExprMat <- lapply(cts, function(ct){
  
  exprMat <- seur[['alra']]@data
  meanExpr <- rowMeans(exprMat[, seur@active.ident == ct])
  
}) %>% do.call("cbind",.) %>% as.data.frame()
colnames(meanExprMat) <- cts

AllMarkers <- read.csv("./result/all.markers.alra.csv",row.names = 1)
DEGs_list <- split(AllMarkers$gene,AllMarkers$cluster)

## load Database ####################

library(OmnipathR)

ppi_form_kinaseextra <- import_kinaseextra_interactions() 
ppi_form_omnipath <- import_omnipath_interactions() 
ppi_form_pathwayextra <- import_pathwayextra_interactions() 
ppi_from_Omnipath <- rbind(ppi_form_kinaseextra,ppi_form_omnipath,ppi_form_pathwayextra) # 101590 pairs
ppiDB <- ppi_from_Omnipath %>% dplyr::filter(n_resources>0) %>% 
  dplyr::select(source_genesymbol, target_genesymbol, n_resources) %>%
  dplyr::rename(source=source_genesymbol, target=target_genesymbol, weight = n_resources) %>%
  dplyr::filter(source %in% rownames(meanExprMat), target %in% rownames(meanExprMat))

LigRecDB <- import_ligrecextra_interactions() # 6456 pair
LigRecDB <- LigRecDB %>% 
  dplyr::select(source_genesymbol, target_genesymbol, n_resources) %>%
  dplyr::rename(source=source_genesymbol, target=target_genesymbol, weight = n_resources) %>%
  dplyr::filter(source %in% unique(c(ppiDB$source,ppiDB$target)), 
                target %in% unique(c(ppiDB$source,ppiDB$target)))

TFTG_from_dorothea <- import_dorothea_interactions() 
TFTG_from_tftarget <- import_tf_target_interactions() 
TFTGDB <- rbind(TFTG_from_dorothea[,-grep("dorothea_level",colnames(TFTG_from_dorothea))],
                TFTG_from_tftarget) # 77182 pairs
TFTGDB <- TFTGDB %>% 
  dplyr::select(source_genesymbol, target_genesymbol, n_resources) %>%
  dplyr::rename(source=source_genesymbol, target=target_genesymbol, weight = n_resources) %>%
  dplyr::filter(source %in% unique(c(ppiDB$source,ppiDB$target)), 
                target %in% unique(c(ppiDB$source,ppiDB$target))) 

Database <- list(LigRecDB = LigRecDB,
                 ppiDB = ppiDB,
                 TFTGDB = TFTGDB)
saveRDS(Database, "./result/database.rds")

## creat undirected weighted graph ##########

V = unique(c(ppiDB$source,ppiDB$target))
adjM <- matrix(0,nrow = length(V),ncol = length(V), dimnames = list(V,V))
for(i in 1:nrow(ppiDB)){
  
  adjM[ppiDB$source[i],ppiDB$target[i]] <- 1/ppiDB$weight[i]
  adjM[ppiDB$target[i],ppiDB$source[i]] <- 1/ppiDB$weight[i]
  
}
adjM <- as(adjM, "dgCMatrix")

G <- igraph::graph_from_adjacency_matrix(adjM, mode='undirected', weighted=TRUE)
# delete loop and multiple edges
G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
# delete isolated nodes
isolates <- which(degree(G, mode = c("all")) == 0) - 1
G <- delete.vertices(G, names(isolates))

Receptors <- LigRecDB$target %>% as.character() %>% unique()
TFs <- TFTGDB$source %>% as.character() %>% unique()
distM <- lapply(Receptors, function(rec){
  
  distances(G, v = rec, to = TFs)
  
}) %>% do.call("rbind",.) %>% as.data.frame()
rownames(distM) <- Receptors
colnames(distM) <- TFs

score.cutoff = 0.2
RecTFDB <- distM
RecTFDB <- cbind(rownames(RecTFDB),RecTFDB)
RecTFDB <- reshape2::melt(RecTFDB, id = "rownames(RecTFDB)")
colnames(RecTFDB) <- c('source','target','score')
RecTFDB <- RecTFDB %>%
  dplyr::filter(score != Inf, score != 0) %>% # 498940 pairs
  dplyr::filter(score <= quantile(score, score.cutoff)) # 48856 pairs

saveRDS(RecTFDB, "./result/RecTFDB.rds")

## main function #######

runMLnet <- function(sender,reciver,rank.cutoff=0.2){
  
  ## get LigRec #################
  
  LigRec <- LigRecDB
  LigRec$LRscore <- apply(LigRec,1,function(x){
    
    meanExprMat[x[1],sender]*meanExprMat[x[2],reciver]
    
  })
  LigRec$rank <- rank(-LigRec$LRscore) 
  LigRec <- LigRec[LigRec$rank<quantile(LigRec$rank,rank.cutoff),]
  LigRec <- LigRec[LigRec$target %in% DEGs_list[[reciver]],]
  
  ## get TFTG #################
  
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
  
  target_gene <- rownames(meanExprMat)
  target_deg <- DEGs_list[[reciver]]
  tryCatch({
    TFTG <- getTFTG(TFTGDB, target_deg, target_gene)
  },
  error=function(e){cat(conditionMessage(e),"\n")})
  tag1 = exists("TFTG")
  if(!tag1) return(NULL)
  
  ## get RecTF ################
  
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
    Recs <- Recs[Recs %in% target_gene]
    
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

## perform ########

df_cellpair <- paste(rep(cts,each=length(cts)), rep(cts,times=length(cts)), sep = "_")
df_cellpair <- lapply(df_cellpair, function(x){strsplit(x,"_")[[1]]}) %>% do.call("rbind",.)
df_cellpair <- df_cellpair[apply(df_cellpair, 1, function(x){x[1]!=x[2]}),]
df_cellpair <- data.frame(df_cellpair)
colnames(df_cellpair) <- c("sender","receiver")

MLnet <- apply(df_cellpair, 1, function(cp){
  
  cat(paste0("sender: ",cp[1],'\n'))
  cat(paste0("receiver: ",cp[2],'\n'))
  res <- runMLnet(cp[1],cp[2])
  cat("###########################\n")
  res
  
})
names(MLnet) <- paste(df_cellpair$sender, df_cellpair$receiver, sep = "_")
str(MLnet)
saveRDS(MLnet, "./result/MLnet.rds")

## plot input #################

rm(list = ls())
MLnet <- readRDS("./result/MLnet.rds")
df_cellpair <- data.frame(sender = lapply(names(MLnet), function(x){strsplit(x,"_")[[1]][1]}) %>% unlist(),
                          receiver = lapply(names(MLnet), function(x){strsplit(x,"_")[[1]][2]}) %>% unlist())

## plot communication network in microenvironment #######

# library
library(Matrix)
library(dplyr)
library(ggplot2)
library(igraph)

# install.packages("plotrix")
library(plotrix)

# function 
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
  # allcolos <- rainbow(length(alltype))
  ChoColorNum <- 1:length(alltype) # sample(1:length(colodb),length(alltype), replace = FALSE)
  
  allcolos <- colodb[ChoColorNum]
  names(allcolos) <- alltype
  
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

# color
colodb <- mycols_seur[c(2,4,3,1)]
show_col(colodb)

# plot input
LRtab <- df_cellpair
colnames(LRtab) <- c('cell_from','cell_to')
LRtab$n <- lapply(MLnet,function(net){
  
  if(!is.null(net)){
    n=nrow(net$LigRec)
  }else{
    n=0
  }
  
}) %>% unlist()

# plot
pdf("./result/cell_comm.pdf",height = 10,width = 20)
DrawCellComm(LRtab,colodb)
dev.off()

## plot Gal-9-Tim-3-x-M2 network ###########

# plot input
lig = "LGALS9"
rec = "HAVCR2"
mlnet_m2 <- list(LigRec = data.frame(source = lig, target = rec),
                 RecTF = mlnet$RecTF[mlnet$RecTF$source %in% rec & mlnet$RecTF$target %in% tfs,],
                 TFTar = mlnet$TFTG[mlnet$TFTG$source %in% tfs & mlnet$TFTG$target %in% m2tgs,])


df_edges <- do.call("rbind",mlnet_m2)
df_nodes <- data.frame(node = c(lig,rec,unique(mlnet_m2$TFTar$source),unique(mlnet_m2$TFTar$target)))
df_nodes$key <- c('Ligand','Receptor',
                  rep('TF',length(intersect(df_nodes$node,tfs))),
                  rep('Target',length(intersect(df_nodes$node,m2tgs))))

# color
library(RColorBrewer)
coul <- brewer.pal(nlevels(as.factor(df_nodes$key)), "Set2")
my_color <- coul[as.numeric(as.factor(df_nodes$key))]

# plot

# pdf("./result/network_reingold.tilford.pdf")

network <- graph_from_edgelist(as.matrix(df_edges[,1:2]))
par(bg="white", mar=c(0,0,0,3))
set.seed(4)
plot(network, 
     vertex.size=12,
     vertex.color=my_color, 
     vertex.label.cex=0.9,
     vertex.label.color="black",
     vertex.frame.color="transparent",
     edge.width=2,           
     edge.arrow.size=0.3,
     edge.arrow.width=0.8,   
     layout = layout.reingold.tilford
)

legend(x=1.1, y=-0.5, 
       legend=levels(as.factor(df_nodes$key)), 
       col = coul , 
       bty = "n", pch=20 , pt.cex = 2, cex = 1,
       text.col="black" , horiz = F)

# dev.off()

#################################
## Enrichment analysis Tim3-up ##
#################################
## library ###########

library(Seurat) 
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# BiocManager::install("ggsci")

library(RColorBrewer)
display.brewer.all()
mycols <- brewer.pal(11, "RdBu")
mycols <- mycols[c(1,4)]

## genelist ############

MLnet <- readRDS("./result/MLnet.rds")
MLnet <- MLnet$Malignant_macrophages
tfs <- MLnet$RecTF$target[MLnet$RecTF$source %in% "HAVCR2"]
tgs <- MLnet$TFTG$target[MLnet$TFTG$source  %in% tfs]
geneList <- tgs

g2s <- toTable(org.Hs.egSYMBOL)
geneList <- g2s$gene_id[match(geneList,g2s$symbol)]
geneList <- na.omit(geneList)

## GO ###########

enrich_gobp <- enrichGO(geneList, 'org.Hs.eg.db', ont="BP", keyType = 'ENTREZID', 
                        minGSSize = 1, pvalueCutoff = 0.99)
str(enrich_gobp) # 6372 genesets

res_enrich_gobp <- enrich_gobp@result
res_enrich_gobp <- res_enrich_gobp[res_enrich_gobp$p.adjust <= 0.05,]

hit_terms <- res_enrich_gobp$Description
hit_terms[grep("macrophage activation",hit_terms,ignore.case = T)] # 3 terms
hit_terms[grep("vessel",hit_terms,ignore.case = T)] # 2 terms
hit_terms[grep("vasculature",hit_terms,ignore.case = T)] # 2 terms

hits <- res_enrich_gobp$Description[grep("macrophage activation|vessel diameter|of vasculature development",
                                         res_enrich_gobp$Description,ignore.case = T)]
res_enrich_gobp_hits <- res_enrich_gobp[res_enrich_gobp$Description %in% hits,]
res_enrich_gobp_hits$GeneRatio <- res_enrich_gobp_hits$Count/343

# pdf("./result/gsea/enrich_gobp.pdf", width = 8, height = 2)
p1 <- ggplot(res_enrich_gobp_hits, aes(x = GeneRatio, y = reorder(Description ,GeneRatio), 
                                 size = Count, fill = p.adjust)) + 
  geom_point(shape = 21) + ylab("GO Term") + theme_bw() +
  scale_fill_continuous(low = mycols[1], high = mycols[2]) + 
  theme(axis.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.box = "horizontal") 
p1
# dev.off()

## KEGG ##############

enrich_kegg <- enrichKEGG(geneList, minGSSize = 1, pvalueCutoff = 0.99)
str(enrich_kegg) # 274 genesets

res_enrich_kegg <- enrich_kegg@result
res_enrich_kegg <- res_enrich_kegg[res_enrich_kegg$p.adjust <= 0.05,]

hit_terms <- res_enrich_kegg$Description
hit_terms[grep("Estrogen|NF-kappa|JAK-STAT",hit_terms,ignore.case = T)] # 3 terms

res_enrich_kegg_hits <- res_enrich_kegg[1:10,]
res_enrich_kegg_hits$GeneRatio <- res_enrich_kegg_hits$Count/316

# pdf("./result/gsea/enrich_kegg.pdf", width = 8, height = 3)
p2 <- ggplot(res_enrich_kegg_hits, aes(x = GeneRatio, y = reorder(Description ,GeneRatio), 
                                 size = Count, fill = p.adjust)) + 
  geom_point(shape = 21) + ylab("KEGG Term") + theme_bw() +
  scale_fill_continuous(low = mycols[1], high = mycols[2]) + 
  theme(axis.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.box = "horizontal") 
p2
# dev.off()

pdf("./result/gsea/enrich.pdf", width = 8, height = 5)
ggpubr::ggarrange(p1,p2,ncol = 1,align = "hv",heights = c(2,3))
dev.off()

###################
## GSEA analysis ##
###################
## library ########

rm(list = ls())

library(Seurat)
library(clusterProfiler)
library(enrichplot)

# BiocManager::install('GSEABase')
library(GSEABase) 

# BiocManager::install("ggsci")
library(ggsci)
mycols <- pal_lancet("lanonc", alpha = 0.7)(9)

## load input ###########

seur <- readRDS("./result/seur.rds")
DefaultAssay(seur) <- 'alra'
degs.all <- FindMarkers(seur, ident.1 = "macrophages", logfc.threshold = 0.01, only.pos = FALSE)

## genelist ############

geneList <-  degs
geneList <- geneList$avg_log2FC
names(geneList) = rownames(degs)
geneList <- geneList[degs$pct.1>=0.1]
geneList <- geneList[order(geneList, decreasing = T)]

## GO ############

gmtfile ='./msigDB/c5.go.bp.v7.4.symbols.gmt'
geneset <- read.gmt(gmtfile)
geneset$term <- geneset$term %>% tolower()
length(unique(geneset$term))
# 7481 genesets

gsea_gobp <- GSEA(geneList, TERM2GENE=geneset, minGSSize = 1, pvalueCutoff = 0.99, verbose=FALSE, seed = 10)
str(gsea_gobp) # 7192 genesets

res_gsea_gobp <- gsea_gobp@result
res_gsea_gobp <- res_gsea_gobp[res_gsea_gobp$p.adjust <= 0.05,]

hit_terms <- res_gsea_gobp$ID
hit_terms[grep("macrophage_activation",hit_terms,ignore.case = T)]
hit_terms[grep("vessel_morphogenesis",hit_terms,ignore.case = T)]
hit_terms[grep("vasculature",hit_terms,ignore.case = T)]

hits <- res_gsea_gobp$ID[grep("macrophage_activation",res_gsea_gobp$ID,ignore.case = T)]
# pdf("./result/gsea/macrophages_activation.pdf", width = 12, height = 8)
gseaplot2(gsea_gobp, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
# dev.off()

hits <- res_gsea_gobp$ID[grep("vessel_morphogenesis|vasculature",res_gsea_gobp$ID,ignore.case = T)]
# pdf("./result/gsea/vessel_development.pdf", width = 12, height = 8)
gseaplot2(gsea_gobp, geneSetID = hits, pvalue_table = T, base_size = 18, color = mycols[1:length(hits)])
# dev.off()
