rm(list=ls())

.libPaths("/home/kzayne/R/x86_64-pc-linux-gnu-library/4.2/")

library(patchwork)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(readxl)
library(HGNChelper)
library(DESeq2)
library(MAST)
library(scCustomize)
library(ggplot2)
library(ggplotify)
library(ggpubr)
library(glmGamPoi)
library(sctransform)
# library(harmony)
library(cowplot)
library(DoubletFinder)
library(Cairo)
library(fgsea)
library(GSA)
library(doBy)
library(EnhancedVolcano)
library(edgeR)
library(RColorBrewer)
# library(infercnv)
# library(monocle3)
# library(metap)
# library(multtest)

setdir <- "/home/kzayne/Salami-Lab-RCC-Organoid-Project/"
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

### Cell markers file
# kid.mrkrs <- read.csv("Updated Kidney Markers original.csv",header = T) #Final kidney markers.csv
# kid.mrkrs <- kid.mrkrs[!kid.mrkrs$cell_name %in% c("Neutrophil","Cancer stem cell"),]
# mrkr.list <- as.list(as.character(kid.mrkrs$Symbol))
# names(mrkr.list) <- kid.mrkrs$cell_name
# 
# for(i in 1:length(mrkr.list)){
#   mrkr.list[[i]] <- unlist(strsplit(mrkr.list[[i]], "," ))
# }

kid.mrkrs <- read.csv("Mohan Normal Kidney Cell Markers.csv",header = T) #Final kidney markers.csv
kid.mrkrs <- kid.mrkrs[!kid.mrkrs$cell_name %in% c("Neutrophil","Cancer stem cell"),]
mrkr.list <- as.list(as.character(kid.mrkrs$Symbol))
names(mrkr.list) <- kid.mrkrs$cell_name

for(i in 1:length(mrkr.list)){
  mrkr.list[[i]] <- unlist(strsplit(mrkr.list[[i]], "," ))
}

# mrkr.list1[['Tumor']] = c(mrkr.list1[['Proximal tubule-1']],'NNMT','CA9','KRT19','KRT18')

# Hallmark cancer
hm.sym <- GSA.read.gmt("/home/kzayne/Salami-Lab-RCC-Organoid-Project/h.all.v7.4.symbols.gmt.txt")    # download all available genesets genesymbols .gmt file from MSigDB
names(hm.sym$genesets) <- hm.sym$geneset.names
hm.sym <- hm.sym$genesets
names(hm.sym) <- gsub('HALLMARK_','',names(hm.sym))


### Load in data
#####
# Read in normal tissue
rcc10n.data <- Read10X(data.dir = '/avatar_data4/11632-TP/10x_analysis_11632-TP/Sample_11632-TP-1/filtered_feature_bc_matrix')
rcc10n <- CreateSeuratObject(counts = rcc10n.data, project = 'RCC10N Tissue')

# Read in tumor 1 tissue
rcc10t1.data <- Read10X(data.dir = '/avatar_data4/11632-TP/10x_analysis_11632-TP/Sample_11632-TP-2/filtered_feature_bc_matrix/')
rcc10t1 <- CreateSeuratObject(counts = rcc10t1.data, project = 'RCC10T1 Tissue')

# Read in tumor 2 tissue
rcc10t2.data <- Read10X(data.dir = '/avatar_data4/11632-TP/10x_analysis_11632-TP/Sample_11632-TP-3/filtered_feature_bc_matrix/')
rcc10t2 <- CreateSeuratObject(counts = rcc10t2.data, project = 'RCC10T2 Tissue')
#####

# Get mitochondrial RNA percentage
rcc10n[['percent.mt']] <- PercentageFeatureSet(rcc10n, pattern = '^MT-')
rcc10t1[['percent.mt']] <- PercentageFeatureSet(rcc10t1, pattern = '^MT-')
rcc10t2[['percent.mt']] <- PercentageFeatureSet(rcc10t2, pattern = '^MT-')

### Generate QC Figures Prior to any filtering
#####
rcc10n.pre <- subset(rcc10n, subset = percent.mt < 10)
rcc10t1.pre <- subset(rcc10t1, subset = percent.mt < 10)
rcc10t2.pre <- subset(rcc10t2, subset = percent.mt < 10)

set.seed(555)
# General workflow for rcc10n
rcc10n.pre <- NormalizeData(rcc10n.pre, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc10n.pre@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# General workflow for rcc10t1
rcc10t1.pre <- NormalizeData(rcc10t1.pre, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc10t1.pre@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# General workflow for rcc10t2
rcc10t2.pre <- NormalizeData(rcc10t2.pre, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc10t2.pre@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# RCC10N doublet finder; 11,318 cells -> 8.8% multiplet rate
rcc10n.pre$multRate <- 0.088 # from 10X based on cells recovered
rcc10n.pre <- NormalizeData(rcc10n.pre, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc10n.pre@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc10n.pre, PCs = 1:20, sct = F)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = F)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)

# ggplot(bcmvn_nsclc, aes(pK, BCmetric, group = 1)) +
#   geom_point() +
#   geom_line()

pK <- bcmvn_nsclc %>% # select pK with max BCmetric
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

# Homotypic doublet proprotion estimate
annotations <- rcc10n.pre@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc10n.pre$multRate)*nrow(rcc10n.pre@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc10n.pre <- doubletFinder(rcc10n.pre,
                             PCs = 1:20,
                             pN = 0.25,
                             pK = pK,
                             nExp = nExp_poi.adj,
                             reuse.pANN = F,
                             sct = F)
colDF <- colnames(rcc10n.pre@meta.data)[grepl('DF.*',colnames(rcc10n.pre@meta.data))]

# Rename doublet identying column
names(rcc10n.pre@meta.data)[names(rcc10n.pre@meta.data) == colDF] <- 'DoubletID'
doubletN <- DimPlot(rcc10n.pre, 
        reduction = 'umap', 
        group.by = 'DoubletID') + ggtitle('RCC10 N')

# RCC10T1 doublet finder; 20,392 cells -> 16% multiplet rate
rcc10t1.pre$multRate <- 0.16 # from 10X based on cells recovered
rcc10t1.pre <- NormalizeData(rcc10t1.pre, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc10t1.pre@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc10t1.pre, PCs = 1:20, sct = F)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = F)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)

# ggplot(bcmvn_nsclc, aes(pK, BCmetric, group = 1)) +
#   geom_point() +
#   geom_line()

pK <- bcmvn_nsclc %>% # select pK with max BCmetric
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

# Homotypic doublet proprotion estimate
annotations <- rcc10t1.pre@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc10t1.pre$multRate)*nrow(rcc10t1.pre@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc10t1.pre <- doubletFinder(rcc10t1.pre,
                        PCs = 1:20,
                        pN = 0.25,
                        pK = pK,
                        nExp = nExp_poi.adj,
                        reuse.pANN = F,
                        sct = F)
colDF <- colnames(rcc10t1.pre@meta.data)[grepl('DF.*',colnames(rcc10t1.pre@meta.data))]

# Rename doublet identying column
names(rcc10t1.pre@meta.data)[names(rcc10t1.pre@meta.data) == colDF] <- 'DoubletID'
doubletT1 <- DimPlot(rcc10t1.pre, 
        reduction = 'umap', 
        group.by = 'DoubletID') + ggtitle('RCC10 T1')

# RCC10T2 doublet finder; 19,929 cells -> 16% multiplet rate
rcc10t2.pre$multRate <- 0.16 # from 10X based on cells recovered
rcc10t2.pre <- NormalizeData(rcc10t2.pre, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc10t2.pre@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc10t2.pre, PCs = 1:20, sct = F)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = F)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)

# ggplot(bcmvn_nsclc, aes(pK, BCmetric, group = 1)) +
#   geom_point() +
#   geom_line()

pK <- bcmvn_nsclc %>% # select pK with max BCmetric
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

# Homotypic doublet proprotion estimate
annotations <- rcc10t2.pre@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc10t2.pre$multRate)*nrow(rcc10t2.pre@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc10t2.pre <- doubletFinder(rcc10t2.pre,
                        PCs = 1:20,
                        pN = 0.25,
                        pK = pK,
                        nExp = nExp_poi.adj,
                        reuse.pANN = F,
                        sct = F)
colDF <- colnames(rcc10t2.pre@meta.data)[grepl('DF.*',colnames(rcc10t2.pre@meta.data))]

# Rename doublet identying column
names(rcc10t2.pre@meta.data)[names(rcc10t2.pre@meta.data) == colDF] <- 'DoubletID'
doubletT2 <- DimPlot(rcc10t2.pre, 
        reduction = 'umap', 
        group.by = 'DoubletID') + ggtitle('RCC10 T2')


### Combine individual objects
rcc10.pre <- merge(rcc10n.pre, y = c(rcc10t1.pre, rcc10t2.pre),
               add.cell.ids = c("rcc10N", "rcc10T1", "rcc10T2")) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc10.pre@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(resolution = 0.8, verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

rcc10.pre$orig.ident <- factor(rcc10.pre$orig.ident,
                           levels = c('RCC10N Tissue','RCC10T1 Tissue','RCC10T2 Tissue'),
                           labels = c('RCC10N','RCC10T1','RCC10T2'))
table(rcc10.pre$orig.ident)
clusterPre <- DimPlot(rcc10.pre, group.by = 'seurat_clusters', split.by = 'orig.ident', label = T)


# Violin plots
violinPre <- VlnPlot(rcc10.pre,
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3,
                     group.by = 'orig.ident')

# Scatter plots
scatterCountMT.pre <- FeatureScatter(rcc10.pre, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident')

# Scatter plots
scatterCountFeature.pre <- FeatureScatter(rcc10.pre, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')

p1o.pre <- DimPlot(object = rcc10.pre, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o.pre <- VlnPlot(object = rcc10.pre, features = "PC_1", group.by = "orig.ident", pt.size = .1)
# p1o | p2o

p3o.pre <- DimPlot(object = rcc10.pre,
               reduction = 'umap',
               pt.size = .1,
               group.by = "orig.ident")
identUMAP <- DimPlot(object = rcc10.pre,
                     reduction = 'umap',
                     pt.size = .1,
                     group.by = 'orig.ident')
doubletUMAP <- DimPlot(object = rcc10.pre,
                       reduction = 'umap',
                       pt.size = .1,
                       group.by = 'DoubletID') + ggtitle('RCC10 Combined')
qcDoublet <- identUMAP / doubletUMAP
# ElbowPlot(rcc10)
qcPlots.pre <- ggarrange(scatterCountMT.pre,scatterCountFeature.pre,p1o.pre,p3o.pre,ncol=2,nrow=2) +
  ggtitle('RCC10 Pre-Filtering QC Plots')

doubletPlots <-ggarrange(doubletN,doubletT1,doubletT2,doubletUMAP,ncol = 2,nrow = 2)

#####

### Perform QC Filtering
#####
# Filter out high MT %
#> RCC10N 11,318 'cells' -> 4047 cells, 330 doublets -> 3717
#> RCC10T1 20,392 'cells' -> 20149 cells, 3035 doublets -> 17,114
#> RCC10T2 19,929 'cells' -> 19756 cells, 2964 doublets -> 16,792
rcc10 <- subset(rcc10.pre,
                     subset = DoubletID == 'Singlet')
#####

### Post QC filtering plots
#####
table(rcc10$orig.ident)
dimPost <- DimPlot(rcc10, group.by = 'seurat_clusters', split.by = 'orig.ident', label = T)


# Violin plots
violinPost <- VlnPlot(rcc10,
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3,
                     group.by = 'orig.ident')

# Scatter plots
scatterCountMT.post <- FeatureScatter(rcc10, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident')

scatterCountFeature.post <- FeatureScatter(rcc10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')

p1o.post <- DimPlot(object = rcc10, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o.post <- VlnPlot(object = rcc10, features = "PC_1", group.by = "orig.ident", pt.size = .1)

p3o.post <- DimPlot(object = rcc10,
               reduction = 'umap',
               pt.size = .1,
               group.by = "orig.ident")
# ElbowPlot(rcc10)
qcPlots.post <- ggarrange(scatterCountMT.post,scatterCountFeature.post,p1o.post,p3o.post,ncol=2,nrow=2)
#####

# Visualization
orig.stim <- DimPlot(rcc10, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim

# NNMT Feature plot
cancer.features <- FeaturePlot(rcc10,
                               features = c('EPOR','CA9'),#c('NDUFA4L2','CA9','VEGFA','EGFR','NNMT','IGF2BP3','CXCL14'),
                               reduction = 'umap',
                               label = T,
                               pt.size = 0.2,
                               repel = T,
                               order = T,
                               min.cutoff = 'q10',
                               max.cutoff = 'q90',
                               split.by = 'orig.ident')
cancer.features

condition <- DimPlot(rcc10, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 3,order=T) + ggtitle("By Tumor Normal")
condition
clusters <- DimPlot(rcc10, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 3,order=T) + ggtitle("By UMAP Clusters")
clusters

clusterOriginPlot <- ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

###############################################################################|
### Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc10) <- "RNA" # SCT
rcc10 <- ScaleData(rcc10,verbose = T,features = rownames(rcc10))
es.max = sctype_score(scRNAseqData = rcc10@assays$RNA$scale.data, #SCT
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc10@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc10@meta.data[rcc10@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc10@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc10@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc10@meta.data$cellassign[rcc10@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
# rcc10$ca9 <- rcc10@assays$RNA$scale.data['CA9',]
# rcc10$cellassign <- ifelse(rcc10$cellassign == 'Proximal tubular cell' & rcc10$ca9 > 0,'Proximal Tubular cell + CA9',rcc10$cellassign)
rcc10@meta.data$cellassign <- ifelse(rcc10@meta.data$cellassign == 'Tumor',
                                    paste0(rcc10@meta.data$cellassign,' ',rcc10@meta.data$seurat_clusters),
                                    rcc10@meta.data$cellassign)

custom1 <- DimPlot(rcc10,
                   reduction = "umap",
                   label = TRUE,
                   repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',
                   pt.size = 0.1,
                   label.size = 5,
                   order = T) +
  ggtitle("Cellassign")
custom1

clusters <- DimPlot(rcc10,
                    reduction = "umap",
                    pt.size = 0.1,
                    split.by='orig.ident',
                    group.by = "seurat_clusters",
                    label = T,
                    repel = T,
                    label.size = 3,
                    order=T) + 
  ggtitle("By UMAP Clusters")

clusterCellassignPlot <- clusters / custom1

custom2 <- DimPlot(rcc10,
                   reduction = "umap",
                   label = TRUE,
                   repel = TRUE,
                   group.by = 'cellassign',
                   pt.size = 0.5,
                   label.size = 3,
                   order = T) +
  ggtitle("Cellassign")
custom2

###############################################################################|
### ssGSEA workflow
# create column for orig.ident and cellassign
rcc10$orig.cellassign <- paste0(rcc10$orig.ident,'.',rcc10$cellassign)

# Get average gene expression by orig cellassign
aggExpr <- AggregateExpression(rcc10,
                               group.by = 'orig.cellassign',
                               normalization.method = 'LogNormalize',
                               return.seurat = T,
                               verbose = F)
#> returns a matrix of logNormalized summed counts by group
# geneCluster <- as.data.frame(geneExpCluster@assays$RNA$scale.data)
# geneCluster <- as.data.frame(geneCluster)


### run escape ssGSEA workflow
require(escape)
require(dittoSeq)
# enrich count data
aggEnrich <- enrichIt(obj = aggExpr@assays$RNA$counts,
                      gene.sets = hm.sym,
                      method = 'ssGSEA',
                      groups = 1000,
                      cores = 2,
                      min.size = 5)
# add enriched counts back to object
aggExpr <- AddMetaData(aggExpr, aggEnrich)
# met.data <- merge(colData(geneExpCluster), aggEnrich, by = "row.names", all=TRUE)
# row.names(met.data) <- met.data$Row.names
# met.data$Row.names <- NULL
# colData(geneExpCluster) <- met.data

# add back annotations of orig.ident and cellasign
aggExpr$orig.ident <- gsub('\\..*','',aggExpr$orig.ident)
aggExpr$cellassign <- gsub('.*\\.','',aggExpr$orig.cellassign)

# heatmap vizualization
colors <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
ssgsea <- dittoHeatmap(aggExpr,
                       genes = NULL,
                       metas = names(aggEnrich),
                       order.by = 'cellassign', #cellassign
                       annot.by = c('cellassign','orig.ident'),
                       cluster_cols = F,
                       cluster_rows = F,
                       heatmap.colors = rev(colors(50)),
                       main = 'RCC10 Pathways by Aggregated Identity and Cellassign')

###############################################################################|
### Agregate by tissue and cluster
rcc10$orig.cluster <- paste0(rcc10$orig.ident,'-',rcc10$seurat_clusters)
rcc10.agg <- AggregateExpression(rcc10,
  return.seurat = T,
  group.by = "orig.cluster",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose = F
)
rcc10.agg <- NormalizeData(rcc10.agg, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc10.pre@var.genes, npcs = 20, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = F)

# Add other columns for break down
rcc10.agg$orig.ident <- gsub('(*.)-.*','\\1',rcc10.agg$orig.cluster)
rcc10.agg$cluster <- gsub('.*-','',rcc10.agg$orig.cluster)


aggDim <- DimPlot(rcc10.agg,
                  reduction = 'umap',
                  label = T,
                  repel = T,
                  pt.size = 0.5,
                  label.size = 3,
                  order = T) +
  NoLegend() +
  ggtitle('RCC10 Aggregated cluster by tissue')

###############################################################################|
### Find Markers
rcc10.join <- JoinLayers(rcc10)
rcc10.allMrkrs <- FindAllMarkers(rcc10.join,
                                min.pct = 0.25,
                                min.diff.pct = 0.25,
                                verbose = F)
top10 <- rcc10.allMrkrs %>% group_by(cluster) %>% top_n(-10, p_val_adj)
# split dataframe into list if you find that convenient
top10.cids <- split(top10$gene, top10$cluster)

### Cell Percentage Plots
pt <- table(rcc10$cellassign, rcc10$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

colourCount = length(unique(rcc10$cellassign))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cellassign <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('rcc10 Percent Cellassign Composition of Samples')

pt2 <- table(rcc10$seurat_clusters, rcc10$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

colourCount = length(unique(rcc10$seurat_clusters))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cluster <- ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('rcc10 Percent Cluster Composition of Samples')

### Generate heatmap of genes
# Get gene list of cell types for heatmap
gene.list <- c()
for(i in 1:nrow(kid.mrkrs)){
  gene.list <- c(gene.list, kid.mrkrs[i,]$Symbol)
}
gene.list <- strsplit(gene.list,',')
gene.list <- unlist(gene.list)  

# Plot heatmap
feature_heatmapn <- DoHeatmap(subset(rcc10,subset = orig.ident == 'RCC10N Tissue'),
                             features = gene.list,#VariableFeatures(rcc9)[1:150],c('CA9','NDUFA4L2','NNMT','VEGFA','HIF1A'),#
                             #cells = 1:500,
                             group.by = 'cellassign',#'seurat_clusters',
                             size = 4,
                             angle = 90) +
  #scale_y_discrete(breaks = c(10,5,5,6,8,8,5,5,6,5,4,4,6,6,10,9,5,8,4,8,5,8,9,4,5,9)) +
  ggtitle('RCC10N Heatmap Tissue')

# Plot heatmap
feature_heatmapt1 <- DoHeatmap(subset(rcc10,subset = orig.ident == 'RCC10T1 Tissue'),
                              features = gene.list,#VariableFeatures(rcc9)[1:150],c('CA9','NDUFA4L2','NNMT','VEGFA','HIF1A'),#
                              #cells = 1:500,
                              group.by = 'cellassign',#'seurat_clusters',
                              size = 4,
                              angle = 90) +
  #scale_y_discrete(breaks = c(10,5,5,6,8,8,5,5,6,5,4,4,6,6,10,9,5,8,4,8,5,8,9,4,5,9)) +
  ggtitle('RCC10T1 Heatmap Tissue')


### Get average gene expression by seurat cluster
geneExpCluster <- AggregateExpression(rcc10,
                                      group.by = 'seurat_clusters',
                                      normalization.method = 'LogNormalize',
                                      return.seurat = T,
                                      verbose = F)
#> returns a matrix of logNormalized summed counts by group
geneCluster <- as.data.frame(geneExpCluster@assays$RNA$scale.data)
geneCluster <- as.data.frame(geneCluster)

### Find gene differences between Normal and Tumor clusters/cell types
# Create column for origin tissue and cluster
rcc10.join$Tumor <- ifelse(rcc10.join$seurat_clusters %in% c(4,9,11),'Tumor 1',
                      ifelse(rcc10.join$seurat_clusters %in% c(7,25,35),'Tumor 2','NonTumor'))
rcc10.join$orig.tumor <- paste0(rcc10.join$orig.ident,'.',rcc10.join$Tumor)
rcc10.join$orig.cluster <- paste0(rcc10$orig.ident,".",rcc10$seurat_clusters)

DimPlot(rcc10.join,
        reduction = 'umap',
        label = T,
        repel = T,
        group.by = 'orig.tumor',
        split.by = 'orig.ident',
        pt.size = 0.1,
        label.size = 3,
        order = T) +
  ggtitle("Identity + tumor group")

# Tumor cluster 0 vs Normal cluster 6
T1.g1_T1.g2 <- FindMarkers(rcc10.join,
                     group.by = 'orig.tumor',
                     ident.1 = 'RCC10T1.Tumor 1',
                     ident.2 = 'RCC10T1.Tumor 2',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
EnhancedVolcano(T1.g1_T1.g2,
                lab = rownames(T1.g1_T1.g2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Tumor Cluster 2 <--Log2FC--> Up in Tumor Cluster 1',
                ylab = 'Adjusted P-value',
                title = 'DEGs in RCC10 T1 Tumor 1 vs Tumor 2',
                pCutoff = 0.05,
                FCcutoff = 1.5)
T1.g1_T1.g2$genes <- rownames(T1.g1_T1.g2)

T1.g1_T1.g2 <- T1.g1_T1.g2 %>% arrange(desc(avg_log2FC))
fold_changes <- T1.g1_T1.g2$avg_log2FC
names(fold_changes) <-T1.g1_T1.g2$genes
T1.g1_T1.g2.gsea <- fgsea(pathways = hm.sym,
                    stats = fold_changes,
                    eps = 0.0,
                    minSize = 15,
                    maxSize = 500)
T1.g1_T1.g2.gsea$comp <- 'T1 Tumor 1 v Tumor 2'

# Tumor cluster 0 vs Normal cluster 6
T2.g1_T2.g2 <- FindMarkers(rcc10.join,
                           group.by = 'orig.tumor',
                           ident.1 = 'RCC10T2.Tumor 1',
                           ident.2 = 'RCC10T2.Tumor 2',
                           min.pct = 0.25,
                           min.diff.pct = 0.25,
                           verbose = F)
EnhancedVolcano(T2.g1_T2.g2,
                lab = rownames(T2.g1_T2.g2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Tumor Cluster 2 <--Log2FC--> Up in Tumor Cluster 1',
                ylab = 'Adjusted P-value',
                title = 'DEGs in RCC10 T2 Tumor 1 vs Tumor 2',
                pCutoff = 0.05,
                FCcutoff = 1.5)
T2.g1_T2.g2$genes <- rownames(T2.g1_T2.g2)

T2.g1_T2.g2 <- T2.g1_T2.g2 %>% arrange(desc(avg_log2FC))
fold_changes <- T2.g1_T2.g2$avg_log2FC
names(fold_changes) <-T2.g1_T2.g2$genes
T2.g1_T2.g2.gsea <- fgsea(pathways = hm.sym,
                          stats = fold_changes,
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)
T2.g1_T2.g2.gsea$comp <- 'T2 Tumor 1 v Tumor 2'

# Tumor cluster 0 vs Normal cluster 6
T1.g1_T2.g1 <- FindMarkers(rcc10.join,
                           group.by = 'orig.tumor',
                           ident.1 = 'RCC10T1.Tumor 1',
                           ident.2 = 'RCC10T2.Tumor 1',
                           min.pct = 0.25,
                           min.diff.pct = 0.25,
                           verbose = F)
EnhancedVolcano(T1.g1_T2.g1,
                lab = rownames(T1.g1_T2.g1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in T2 <--Log2FC--> Up in T1',
                ylab = 'Adjusted P-value',
                title = 'DEGs in RCC10 T1 Tumor 1 vs T2 Tumor 1',
                pCutoff = 0.05,
                FCcutoff = 1.5)
T1.g1_T2.g1$genes <- rownames(T1.g1_T2.g1)

T1.g1_T2.g1 <- T1.g1_T2.g1 %>% arrange(desc(avg_log2FC))
fold_changes <- T1.g1_T2.g1$avg_log2FC
names(fold_changes) <-T1.g1_T2.g1$genes
T1.g1_T2.g1.gsea <- fgsea(pathways = hm.sym,
                          stats = fold_changes,
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)
T1.g1_T2.g1.gsea$comp <- 'T1 Tumor 1 v T2 Tumor 1'

# Tumor cluster 0 vs Normal cluster 6
T1.g2_T2.g2 <- FindMarkers(rcc10.join,
                           group.by = 'orig.tumor',
                           ident.1 = 'RCC10T1.Tumor 2',
                           ident.2 = 'RCC10T2.Tumor 2',
                           min.pct = 0.25,
                           min.diff.pct = 0.25,
                           verbose = F)
EnhancedVolcano(T1.g2_T2.g2,
                lab = rownames(T1.g2_T2.g2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in T2 <--Log2FC--> Up in T1',
                ylab = 'Adjusted P-value',
                title = 'DEGs in RCC10 T1 Tumor 2 vs T2 Tumor 2',
                pCutoff = 0.05,
                FCcutoff = 1.5)
T1.g2_T2.g2$genes <- rownames(T1.g2_T2.g2)

T1.g2_T2.g2 <- T1.g2_T2.g2 %>% arrange(desc(avg_log2FC))
fold_changes <- T1.g2_T2.g2$avg_log2FC
names(fold_changes) <-T1.g2_T2.g2$genes
T1.g2_T2.g2.gsea <- fgsea(pathways = hm.sym,
                          stats = fold_changes,
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)
T1.g2_T2.g2.gsea$comp <- 'T1 Tumor 2 v T2 Tumor 2'

###############################################################################|
### Combined FGSEA bubble plot
#####
tops <- rbind(T1.g1_T1.g2.gsea,T2.g1_T2.g2.gsea, T1.g1_T2.g1.gsea, T1.g2_T2.g2.gsea)
tops <- tops[tops$padj < 0.05,]
tops <- tops[,c(1:7,9)]
tops <- tops[!grepl("HALLMARK",tops$pathway),]
tops$pathway <- gsub("HALLMARK_","",tops$pathway)
tops$pathway <- gsub("_"," ",tops$pathway)
#tops <- tops[!grepl("ESTROGEN",tops$pathway),]
tops <- tops[order(tops$padj,decreasing = F),]
tops$pathway <- factor(tops$pathway,levels = unique(c(tops$pathway)))

tops$comp <- factor(tops$comp)

tops$NES <- as.numeric(round(tops$NES,digits=2))
tops$Direction <- factor(ifelse(tops$padj < 0.05 & tops$NES>0,"Enriched in Tumor 1",
                                ifelse(tops$padj < 0.05 & tops$NES<0,"Enriched in Tumor 2","Not significant")),
                         levels = c("Enriched in Tumor 1", "Enriched in Tumor 2","Not significant"))
tops <- tops[!tops$Direction=="Not significant",]
tops$absNES <- abs(tops$NES)

# ggplot Volcano plots plus hallmark plot
fill <- c("red","#3182bd")

g1 <- ggplot(tops, aes(x = comp, y = pathway, size = absNES, fill = Direction)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_manual(values =fill) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text=element_text(size=10,colour = "black")) + scale_y_discrete(limits = rev(levels(factor(tops$pathway)))) +
  scale_size(range = c(3,8),breaks = c(1,1.5,2)) +
  
  # ggtitle("edgeR common DEGs") +
  xlab("")+ylab("")+theme_bw() +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 45, vjust = 0, hjust = 0,face="bold",size=10,colour = "black")) +
  theme(axis.text.y = element_text(size =  10,colour = "black")) +
  theme(axis.text=element_text(size=20), legend.text=element_text(size=12))+
  guides(size=guide_legend(title="Normalized Enrichment Score (NES)", title.theme = element_text(
    size = 12,
    colour = "black",
    face = "bold",
    angle = 0))) +
  guides(fill=guide_legend(title="Direction", override.aes = list(size=10), title.theme = element_text(
    size = 12,
    colour = "black",
    face = "bold",
    angle = 0))) #+ coord_equal(2/6)
g1
