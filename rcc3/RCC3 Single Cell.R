rm(list=ls())

myPaths <- .libPaths()
myPaths <- c(myPaths,"/home/kzayne/R/x86_64-pc-linux-gnu-library/R-4.2/")
myPaths <- c(myPaths[length(myPaths)],myPaths[1:length(myPaths)-1])
.libPaths(myPaths)
#.libPaths("/home/kzayne/R/x86_64-pc-linux-gnu-library/R-4.2/")

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

#####
### Load in data
# pathways: .67 = '/mnt/DATA2/ccRCC... ; .80 = /avatar_Data2/ccRCC...
# Read in normal tissue
rcc3n.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/RCC3/SI_32604/filtered_feature_bc_matrix")
rcc3n <- CreateSeuratObject(counts = rcc3n.data,project = "RCC3N")

# Read in tumor 1 tissue
rcc3t1.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/RCC3/SI_32602/filtered_feature_bc_matrix")
rcc3t1 <- CreateSeuratObject(counts = rcc3t1.data,project = "RCC3T1")

# Read in tumor 2 tissue
rcc3t2.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/RCC3/SI_32603/filtered_feature_bc_matrix")
rcc3t2 <- CreateSeuratObject(counts = rcc3t2.data,project = "RCC3T2")

# Read in pre-treat organoid
rcc3t.org.pre.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC3T_org_nerat_treat/outs/filtered_feature_bc_matrix')
rcc3t.org.pre <- CreateSeuratObject(counts = rcc3t.org.pre.data, project = 'RCC3T Org PreTreat')

# Read in post-treat organoid
rcc3t.org.post.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC3T_org_nerat_treat_no/outs/filtered_feature_bc_matrix')
rcc3t.org.post <- CreateSeuratObject(counts = rcc3t.org.post.data, project = 'RCC3T Org PostTreat')

### Perform Mito Filtering
# Get mitochondrial RNA percentage
rcc3n[['percent.mt']] <- PercentageFeatureSet(rcc3n, pattern = '^MT-')
rcc3t1[['percent.mt']] <- PercentageFeatureSet(rcc3t1, pattern = '^MT-')
rcc3t2[['percent.mt']] <- PercentageFeatureSet(rcc3t2, pattern = '^MT-')
rcc3t.org.pre[['percent.mt']] <- PercentageFeatureSet(rcc3t.org.pre, pattern = '^MT-')
rcc3t.org.post[['percent.mt']] <- PercentageFeatureSet(rcc3t.org.post, pattern = '^MT-')

# Filter out high MT %
#> RCC3N 3272 'cells' -> 1013 cells
#> RCC3T1 2470 'cells' -> 1688 cells
#> RCC3T2 3367 'cells' -> 2722 cells
#> RCC3T Org Pre 1416 'cells' -> 1376 cells
#> RCC3T Org Post 303 'cells' -> 270
rcc3n.filt <- subset(rcc3n,
                     subset = #nCount_RNA > 800 &
                       #nFeature_RNA > 500 &
                       percent.mt < 10)
rcc3t1.filt <- subset(rcc3t1,
                      subset = percent.mt < 10)
rcc3t2.filt <- subset(rcc3t2,
                      subset = percent.mt < 10)
rcc3t.org.pre.filt <- subset(rcc3t.org.pre,
                             subset = percent.mt < 10)
rcc3t.org.post.filt <- subset(rcc3t.org.post,
                             subset = percent.mt < 10)


###############################################################################|
### Run Doublet Finder
#####
# RCC3N doublet finder; 1013 cells -> 1006 cells
rcc3n.filt$multRate <- 0.008 # from 10X based on cells recovered
rcc3n.filt <- NormalizeData(rcc3n.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc3n@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc3n.filt, PCs = 1:20, sct = F)
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
annotations <- rcc3n.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc3n.filt$multRate)*nrow(rcc3n.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc3n.filt <- doubletFinder(rcc3n.filt,
                            PCs = 1:20,
                            pN = 0.25,
                            pK = pK,
                            nExp = nExp_poi.adj,
                            reuse.pANN = F,
                            sct = F)
colDF <- colnames(rcc3n.filt@meta.data)[grepl('DF.*',colnames(rcc3n.filt@meta.data))]

# Rename doublet identying column
names(rcc3n.filt@meta.data)[names(rcc3n.filt@meta.data) == colDF] <- 'DoubletID'
DimPlot(rcc3n.filt, 
        reduction = 'umap', 
        group.by = 'DoubletID')

# Filter out doublets
rcc3n.filt <- subset(rcc3n.filt, subset = DoubletID == 'Singlet')

# RCC3T1 Doublet Finder; 1688 cells -> 1677 cells
rcc3t1.filt$multRate <- 0.008 # from 10X based on cells recovered
rcc3t1.filt <- NormalizeData(rcc3t1.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc3t1@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc3t1.filt, PCs = 1:20, sct = F)
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
annotations <- rcc3t1.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc3t1.filt$multRate)*nrow(rcc3t1.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc3t1.filt <- doubletFinder(rcc3t1.filt,
                             PCs = 1:20,
                             pN = 0.25,
                             pK = pK,
                             nExp = nExp_poi.adj,
                             reuse.pANN = F,
                             sct = F)
colDF <- colnames(rcc3t1.filt@meta.data)[grepl('DF.*',colnames(rcc3t1.filt@meta.data))]

# Rename doublet identying column
names(rcc3t1.filt@meta.data)[names(rcc3t1.filt@meta.data) == colDF] <- 'DoubletID'
# DimPlot(rcc3t1.filt, 
#         reduction = 'umap', 
#         group.by = 'DoubletID')

# Filter out doublets
rcc3t1.filt <- subset(rcc3t1.filt, subset = DoubletID == 'Singlet')

# RCC3T2 Doublet Finder; 2722 cells -> 2686 cells
rcc3t2.filt$multRate <- 0.016 # from 10X based on cells recovered
rcc3t2.filt <- NormalizeData(rcc3t2.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc3t2@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc3t2.filt, PCs = 1:20, sct = F)
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
annotations <- rcc3t2.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc3t2.filt$multRate)*nrow(rcc3t2.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc3t2.filt <- doubletFinder(rcc3t2.filt,
                             PCs = 1:20,
                             pN = 0.25,
                             pK = pK,
                             nExp = nExp_poi.adj,
                             reuse.pANN = F,
                             sct = F)
colDF <- colnames(rcc3t2.filt@meta.data)[grepl('DF.*',colnames(rcc3t2.filt@meta.data))]

# Rename doublet identying column
names(rcc3t2.filt@meta.data)[names(rcc3t2.filt@meta.data) == colDF] <- 'DoubletID'
# DimPlot(rcc3t2.filt, 
#         reduction = 'umap', 
#         group.by = 'DoubletID')

# Filter out doublets
rcc3t2.filt <- subset(rcc3t2.filt, subset = DoubletID == 'Singlet')

# RCC3T Pre Treat Organoid Doublet Finder; 1376 cells -> 1369 cells
rcc3t.org.pre.filt$multRate <- 0.008 # from 10X based on cells recovered
rcc3t.org.pre.filt <- NormalizeData(rcc3t.org.pre.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc3t2@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc3t.org.pre.filt, PCs = 1:20, sct = F)
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
annotations <- rcc3t.org.pre.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc3t.org.pre.filt$multRate)*nrow(rcc3t.org.pre.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc3t.org.pre.filt <- doubletFinder(rcc3t.org.pre.filt,
                             PCs = 1:20,
                             pN = 0.25,
                             pK = pK,
                             nExp = nExp_poi.adj,
                             reuse.pANN = F,
                             sct = F)
colDF <- colnames(rcc3t.org.pre.filt@meta.data)[grepl('DF.*',colnames(rcc3t.org.pre.filt@meta.data))]

# Rename doublet identying column
names(rcc3t.org.pre.filt@meta.data)[names(rcc3t.org.pre.filt@meta.data) == colDF] <- 'DoubletID'
# DimPlot(rcc3t.org.pre.filt, 
#         reduction = 'umap', 
#         group.by = 'DoubletID')

# Filter out doublets
rcc3t.org.pre.filt <- subset(rcc3t.org.pre.filt, subset = DoubletID == 'Singlet')

# RCC3T Post Treat Organoid Doublet Finder; 270 cells -> 269 cells
rcc3t.org.post.filt$multRate <- 0.004 # from 10X based on cells recovered
rcc3t.org.post.filt <- NormalizeData(rcc3t.org.post.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc3t2@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc3t.org.post.filt, PCs = 1:20, sct = F)
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
annotations <- rcc3t.org.post.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc3t.org.post.filt$multRate)*nrow(rcc3t.org.post.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc3t.org.post.filt <- doubletFinder(rcc3t.org.post.filt,
                             PCs = 1:20,
                             pN = 0.25,
                             pK = pK,
                             nExp = nExp_poi.adj,
                             reuse.pANN = F,
                             sct = F)
colDF <- colnames(rcc3t.org.post.filt@meta.data)[grepl('DF.*',colnames(rcc3t.org.post.filt@meta.data))]

# Rename doublet identying column
names(rcc3t.org.post.filt@meta.data)[names(rcc3t.org.post.filt@meta.data) == colDF] <- 'DoubletID'
# DimPlot(rcc3t.org.post.filt, 
#         reduction = 'umap', 
#         group.by = 'DoubletID')

# Filter out doublets
rcc3t.org.post.filt <- subset(rcc3t.org.post.filt, subset = DoubletID == 'Singlet')
#####

###############################################################################|
### Combine all RCC3 objects
# Merger normal and tissue data
set.seed(555)
rcc3 <- merge(rcc3n.filt, y = c(rcc3t1.filt,rcc3t2.filt,rcc3t.org.pre.filt,rcc3t.org.post.filt),
              add.cell.ids = c("RCC3N", "RCC3T1","RCC3T2",'RCC3T Org PreTreat','RCC3T Org PostTreat')) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc3@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(resolution = 0.8, verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

rcc3$orig.ident <- factor(rcc3$orig.ident,
                          levels = c('RCC3N','RCC3T1','RCC3T2','RCC3T Org PreTreat','RCC3T Org PostTreat'),
                          labels = c('RCC3N','RCC3T1','RCC3T2','RCC3T Org PreTreat','RCC3T Org PostTreat'))
table(rcc3$orig.ident)
dimPre <- DimPlot(rcc3, group.by = 'seurat_clusters', split.by = 'orig.ident', label = T)


# Violin plots
violinPre <- VlnPlot(rcc3,
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3,
                     group.by = 'orig.ident')

# Scatter plots
plot1 <- FeatureScatter(rcc3, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident')
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')
plot2

p1o <- DimPlot(object = rcc3, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc3, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

p3o <- DimPlot(object = rcc3,
               reduction = 'umap',
               pt.size = .1,
               group.by = "orig.ident")
# ElbowPlot(rcc3)
qcPlots <- ggarrange(plot1,plot2,p1o,p3o,ncol=2,nrow=2)

# Visualization
orig.stim <- DimPlot(rcc3, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim


condition <- DimPlot(rcc3, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 3,order=T) + ggtitle("By Tumor Normal")
condition
clusters <- DimPlot(rcc3, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 3,order=T) + ggtitle("By UMAP Clusters")
clusters

clusterOriginPlot <- ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

###############################################################################|
### Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc3) <- "RNA" # SCT
rcc3 <- ScaleData(rcc3,verbose = T,features = rownames(rcc3))
es.max = sctype_score(scRNAseqData = rcc3@assays$RNA$scale.data, #SCT
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc3@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc3@meta.data[rcc3@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc3@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc3@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc3@meta.data$cellassign[rcc3@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
# rcc3$ca9 <- rcc3@assays$RNA$scale.data['CA9',]
# rcc3$cellassign <- ifelse(rcc3$cellassign == 'Proximal tubular cell' & rcc3$ca9 > 0,'Proximal Tubular cell + CA9',rcc3$cellassign)

# Get gene list of cell types for heatmap
gene.list <- c()
for(i in 1:nrow(kid.mrkrs)){
  gene.list <- c(gene.list, kid.mrkrs[i,]$Symbol)
}
gene.list <- strsplit(gene.list,',')
gene.list <- unlist(gene.list)  

### Feature plot - swap out genes as needed
# Change labels for cellassign
Idents(rcc3) <- 'cellassign'

# Plot features - Tissue + organoids
cancer.features <- FeaturePlot(rcc3,
                               features = 'IGF2BP3',#c('NDUFA4L2','CA9','VEGFA','EGFR','NNMT','IGF2BP3'),
                               reduction = 'umap',
                               label = T,
                               repel = T,
                               order = T,
                               min.cutoff = 'q10',
                               max.cutoff = 'q90',
                               split.by = 'orig.ident',
                               cols = c('lightgray','red'))
cancer.features

# Plot features - Tissue only
cancer.features <- FeaturePlot(subset(rcc3,subset = orig.ident %in% c('RCC3N','RCC3T1','RCC3T2')),
                               features = 'NDUFA4L2',#c('NDUFA4L2','CA9','VEGFA','EGFR','NNMT','IGF2BP3'),
                               reduction = 'umap',
                               label = T,
                               repel = T,
                               order = T,
                               min.cutoff = 'q10',
                               max.cutoff = 'q90',
                               split.by = 'orig.ident',
                               cols = c('lightgray','red'))
cancer.features
# Reset idents to clusters
Idents(rcc3) <- 'seurat_clusters'

### Plot heatmap
feature_heatmap <- DoHeatmap(rcc3,
                             features = gene.list,#VariableFeatures(rcc3)[1:150],c('CA9','NDUFA4L2','NNMT','VEGFA','HIF1A'),#
                             #cells = 1:500,
                             group.by = 'cellassign',#'seurat_clusters',
                             size = 4,
                             angle = 90) +
  #scale_y_discrete(breaks = c(10,5,5,6,8,8,5,5,6,5,4,4,6,6,10,9,5,8,4,8,5,8,9,4,5,9)) +
  ggtitle('RCC3 Heatmap Tissue')

custom1 <- DimPlot(rcc3,
                   reduction = "umap",
                   label = TRUE,
                   repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',
                   pt.size = 0.5,
                   label.size = 3,
                   order = T) +
  ggtitle("Cellassign")
custom1

clusters <- DimPlot(rcc3,
                    reduction = "umap",
                    pt.size = 0.5,
                    split.by='orig.ident',
                    group.by = "seurat_clusters",
                    label = T,
                    repel = T,
                    label.size = 3,
                    order=T) + 
  ggtitle("By UMAP Clusters")
clusters

clusterCellassignPlot <- clusters / custom1

custom2 <- DimPlot(rcc3,
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
### Find Markers
rcc3.join <- JoinLayers(rcc3)
rcc3.allMrkrs <- FindAllMarkers(rcc3.join,
                                min.pct = 0.25,
                                min.diff.pct = 0.25,
                                verbose = F)
top10 <- rcc3.allMrkrs %>% group_by(cluster) %>% top_n(-10, p_val_adj)
# split dataframe into list if you find that convenient
top10.cids <- split(top10$gene, top10$cluster)

### Cell Percentage Plots
rcc3@meta.data$cellassign <- ifelse(rcc3@meta.data$cellassign == 'Tumor',
                                    paste0(rcc3@meta.data$cellassign,' ',rcc3@meta.data$seurat_clusters),
                                    rcc3@meta.data$cellassign)
pt <- table(rcc3$cellassign, rcc3$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

colourCount = length(unique(rcc3$cellassign))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cellassign <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('RCC3 Percent Cellassign Composition of Samples')

pt2 <- table(rcc3$seurat_clusters, rcc3$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

colourCount = length(unique(rcc3$seurat_clusters))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cluster <- ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('RCC3 Percent Cluster Composition of Samples')

### Find gene differences between Normal and Tumor clusters/cell types
# Tumor 1 vs Normal
t1_n <- FindMarkers(rcc3.join,
                    group.by = 'orig.ident',
                    ident.1 = 'RCC3T1',
                    ident.2 = 'RCC3N',
                    min.pct = 0.25,
                    min.diff.pct = 0.25,
                    verbose = F)
EnhancedVolcano(t1_n,
                lab = rownames(t1_n),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Normal <--Log2FC--> Up in Tumor 1',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor 1 cells vs Normal cells',
                pCutoff = 0.05,
                FCcutoff = 1.5)

t1_n$genes <- rownames(t1_n)
t1_n.mrkrs <- t1_n %>% arrange(desc(avg_log2FC))
fold_changes <- t1_n.mrkrs$avg_log2FC
names(fold_changes) <-t1_n.mrkrs$genes
t1_n.gsea <- fgsea(pathways = hm.sym,
                   stats = fold_changes,
                   eps = 0.0,
                   minSize = 15,
                   maxSize = 500)
t1_n.gsea$comp <- 'Tumor 1 v Normal'

# Tumor 2 vs Normal
t2_n <- FindMarkers(rcc3.join,
                    group.by = 'orig.ident',
                    ident.1 = 'RCC3T2',
                    ident.2 = 'RCC3N',
                    min.pct = 0.25,
                    min.diff.pct = 0.25,
                    verbose = F)
EnhancedVolcano(t2_n,
                lab = rownames(t2_n),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Normal <--Log2FC--> Up in Tumor 2',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor 2 cells vs Normal cells',
                pCutoff = 0.05,
                FCcutoff = 1.5)
t2_n$genes <- rownames(t2_n)

t2_n.mrkrs <- t2_n %>% arrange(desc(avg_log2FC))
fold_changes <- t2_n.mrkrs$avg_log2FC
names(fold_changes) <-t2_n.mrkrs$genes
t2_n.gsea <- fgsea(pathways = hm.sym,
                   stats = fold_changes,
                   eps = 0.0,
                   minSize = 15,
                   maxSize = 500)
t2_n.gsea$comp <- 'Tumor 2 v Normal'

# Tumor 1 vs Normal - This run finds no DEGs
t1_t2 <- FindMarkers(rcc3.join,
                     group.by = 'orig.ident',
                     ident.1 = 'RCC3T1',
                     ident.2 = 'RCC3T2',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
EnhancedVolcano(t1_t2,
                lab = rownames(t1_t2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Tumor 2 <--Log2FC--> Up in Tumor 1',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor 1 cells vs Tumor 2 cells',
                pCutoff = 0.05,
                FCcutoff = 1.5)

t1_t2$genes <- rownames(t1_t2)
t1_t2.mrkrs <- t1_t2 %>% arrange(desc(avg_log2FC))
fold_changes <- t1_t2.mrkrs$avg_log2FC
names(fold_changes) <-t1_t2.mrkrs$genes
t1_t2.gsea <- fgsea(pathways = hm.sym,
                    stats = fold_changes,
                    eps = 0.0,
                    minSize = 15,
                    maxSize = 500)
t1_t2.gsea$comp <- 'Tumor 1 v Tumor 2'

### Investigate differences between organoids and normal tissue
# PreTreat vs Normal
preO_n <- FindMarkers(rcc3.join,
                    group.by = 'orig.ident',
                    ident.1 = 'RCC3T Org PreTreat',
                    ident.2 = 'RCC3N',
                    min.pct = 0.25,
                    min.diff.pct = 0.25,
                    verbose = F)
EnhancedVolcano(preO_n,
                lab = rownames(preO_n),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Normal <--Log2FC--> Up in PreTreat',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Pre Treat Organoid cells vs Normal tissue cells',
                pCutoff = 0.05,
                FCcutoff = 1.5)

preO_n$genes <- rownames(preO_n)
preO_n.mrkrs <- preO_n %>% arrange(desc(avg_log2FC))
fold_changes <- preO_n.mrkrs$avg_log2FC
names(fold_changes) <-preO_n.mrkrs$genes
preO_n.gsea <- fgsea(pathways = hm.sym,
                   stats = fold_changes,
                   eps = 0.0,
                   minSize = 15,
                   maxSize = 500)
preO_n.gsea$comp <- 'PreTreat Org v Normal'

# PostTreat organoid vs Normal
postO_n <- FindMarkers(rcc3.join,
                    group.by = 'orig.ident',
                    ident.1 = 'RCC3T Org PostTreat',
                    ident.2 = 'RCC3N',
                    min.pct = 0.25,
                    min.diff.pct = 0.25,
                    verbose = F)
EnhancedVolcano(postO_n,
                lab = rownames(postO_n),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Normal <--Log2FC--> Up in Post Treat',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Post Treat Organoid cells vs Normal cells',
                pCutoff = 0.05,
                FCcutoff = 1.5)
postO_n$genes <- rownames(postO_n)

postO_n.mrkrs <- postO_n %>% arrange(desc(avg_log2FC))
fold_changes <- postO_n.mrkrs$avg_log2FC
names(fold_changes) <-postO_n.mrkrs$genes
postO_n.gsea <- fgsea(pathways = hm.sym,
                   stats = fold_changes,
                   eps = 0.0,
                   minSize = 15,
                   maxSize = 500)
postO_n.gsea$comp <- 'PostTreat Org v Normal'

# PostTreat vs PreTreat
preO_postO <- FindMarkers(rcc3.join,
                     group.by = 'orig.ident',
                     ident.1 = 'RCC3T Org PostTreat',
                     ident.2 = 'RCC3T Org PreTreat',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
EnhancedVolcano(preO_postO,
                lab = rownames(preO_postO),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in PreTreat <--Log2FC--> Up in PostTreat',
                ylab = 'Adjusted P-value',
                title = 'DEGs in PostTreat Organoid cells vs PreTreat Organoid cells',
                pCutoff = 0.05,
                FCcutoff = 1.5)

preO_postO$genes <- rownames(preO_postO)
preO_postO.mrkrs <- preO_postO %>% arrange(desc(avg_log2FC))
fold_changes <- preO_postO.mrkrs$avg_log2FC
names(fold_changes) <-preO_postO.mrkrs$genes
preO_postO.gsea <- fgsea(pathways = hm.sym,
                    stats = fold_changes,
                    eps = 0.0,
                    minSize = 15,
                    maxSize = 500)
preO_postO.gsea$comp <- 'PostTreat Org v PreTreat Org'


### Investigate tumor cluster differences in T1 and T2
# Create column for origin tissue and cluster
rcc3$orig.cluster <- paste0(rcc3$orig.ident,".",rcc3$seurat_clusters)
rcc3$orig.cellassign <- paste0(rcc3$orig.ident,'.',rcc3$cellassign)

# Get average gene expression by orig cellassign
aggExpr <- AggregateExpression(rcc3,
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
dittoHeatmap(aggExpr,
             genes = NULL,
             metas = names(aggEnrich),
             order.by = 'orig.ident',
             annot.by = c('cellassign','orig.ident'),
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = rev(colors(50)),
             main = 'RCC3 Pathways by Aggregated Identity and Cellassign')


# Tumor 1 Cluster 13 vs Cluster 10
t1.c13_c10 <- FindMarkers(rcc1.join,
                          group.by = 'orig.cluster',
                          ident.1 = 'RCC1T1.13',
                          ident.2 = 'RCC1T1.10',
                          min.pct = 0.25,
                          min.diff.pct = 0.25,
                          verbose = F)
EnhancedVolcano(t1.c13_c10,
                lab = rownames(t1.c13_c10),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Cluster10 <--Log2FC--> Up in Cluster13',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor 1 Tumor Cluster 13 vs Tumor Cluster 10',
                pCutoff = 0.05,
                FCcutoff = 1.5)

t1.c13_c10$genes <- rownames(t1.c13_c10)
t1.c13_c10.mrkrs <- t1.c13_c10 %>% arrange(desc(avg_log2FC))
fold_changes <- t1.c13_c10.mrkrs$avg_log2FC
names(fold_changes) <-t1.c13_c10.mrkrs$genes
t1.c13_c10.gsea <- fgsea(pathways = hm.sym,
                         stats = fold_changes,
                         eps = 0.0,
                         minSize = 15,
                         maxSize = 500)
t1.c13_c10.gsea$comp <- 'T1 Cluster 13 v Cluster 10'

# Tumor 2 Cluster 13 vs Cluster 10
t2.c13_c10 <- FindMarkers(rcc1.join,
                          group.by = 'orig.cluster',
                          ident.1 = 'RCC1T2.13',
                          ident.2 = 'RCC1T2.10',
                          min.pct = 0.25,
                          min.diff.pct = 0.25,
                          verbose = F)
EnhancedVolcano(t2.c13_c10,
                lab = rownames(t2.c13_c10),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Cluster10 <--Log2FC--> Up in Cluster13',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor 2 Tumor Cluster 13 vs Tumor Cluster 10',
                pCutoff = 0.05,
                FCcutoff = 1.5)

t2.c13_c10$genes <- rownames(t2.c13_c10)
t2.c13_c10.mrkrs <- t2.c13_c10 %>% arrange(desc(avg_log2FC))
fold_changes <- t2.c13_c10.mrkrs$avg_log2FC
names(fold_changes) <-t2.c13_c10.mrkrs$genes
t2.c13_c10.gsea <- fgsea(pathways = hm.sym,
                         stats = fold_changes,
                         eps = 0.0,
                         minSize = 15,
                         maxSize = 500)
t2.c13_c10.gsea$comp <- 'T2 Cluster 13 v Cluster 10'

###############################################################################|
### Combined FGSEA bubble plot
#####
tops <- rbind(t1.c13_c10.gsea,t2.c13_c10.gsea)
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
tops$Direction <- factor(ifelse(tops$padj < 0.05 & tops$NES>0,"Enriched in Cluster 13",
                                ifelse(tops$padj < 0.05 & tops$NES<0,"Enriched in Cluster 10","Not significant")),
                         levels = c("Enriched in Cluster 13", "Enriched in Cluster 10","Not significant"))
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

#####