rm(list=ls())

.libPaths("/home/kzayne/R/x86_64-pc-linux-gnu-library/R-4.2/")

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
library(harmony)
library(cowplot)
library(DoubletFinder)
library(Cairo)
library(fgsea)
library(GSA)
library(doBy)
library(EnhancedVolcano)
library(edgeR)
library(RColorBrewer)
library(monocle3)
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
                               features = c('CXCL14','CA9'),#c('NDUFA4L2','CA9','VEGFA','EGFR','NNMT','IGF2BP3','CXCL14'),
                               reduction = 'umap',
                               label = T,
                               pt.size = 0.2,
                               repel = T,
                               order = T,
                               min.cutoff = 1.5,
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
rcc10@meta.data$cellassign <- ifelse(rcc10@meta.data$cellassign == 'Tumor',
                                    paste0(rcc10@meta.data$cellassign,' ',rcc10@meta.data$seurat_clusters),
                                    rcc10@meta.data$cellassign)
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
rcc10.allMrkrs %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(rcc10,group.by = 'cellassign', features = top10$gene) + NoLegend()

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
rcc10.join$orig.cluster <- paste0(rcc10$orig.ident,".",rcc10$seurat_clusters)
# Tumor cluster 0 vs Normal cluster 6
t0_n6 <- FindMarkers(rcc10.join,
                     group.by = 'orig.cluster',
                     ident.1 = 'rcc10T.0',
                     ident.2 = 'rcc10N.6',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
EnhancedVolcano(t0_n6,
                lab = rownames(t0_n6),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Normal Prox-Tub <--Log2FC--> Up in Tumor 0',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor cells in Tumor Cluster 0 vs Normal Proximal Tubule 2',
                pCutoff = 0.05,
                FCcutoff = 1.5)