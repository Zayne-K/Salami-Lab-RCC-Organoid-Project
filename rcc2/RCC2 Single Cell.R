rm(list=ls())

myPaths <- .libPaths()
myPaths <- c(myPaths,"/home/kzayne/R/x86_64-pc-linux-gnu-library/R-4.2/")
myPaths <- c(myPaths[length(myPaths)],myPaths[1:length(myPaths)-1])
.libPaths(myPaths)

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
rcc2n.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/10x_analysis_5964-SS/Sample_5964-SS-1_ACGACTACCA-TTAGGGTCGT/filtered_feature_bc_matrix")
rcc2n <- CreateSeuratObject(counts = rcc2n.data,project = "RCC2N")

# Read in tumor 1 tissue
rcc2t1.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/10x_analysis_5964-SS/Sample_5964-SS-2_CCCAGCTTCT-GTTTGGTGTC/filtered_feature_bc_matrix")
rcc2t1 <- CreateSeuratObject(counts = rcc2t1.data,project = "RCC2T1")

# Read in tumor 2 tissue
rcc2t2.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/10x_analysis_5964-SS/Sample_5964-SS-3_TTGAGAGTCA-CTACCAGGTT/filtered_feature_bc_matrix")
rcc2t2 <- CreateSeuratObject(counts = rcc2t2.data,project = "RCC2T2")

### Perform Mito Filtering
# Get mitochondrial RNA percentage
rcc2n[['percent.mt']] <- PercentageFeatureSet(rcc2n, pattern = '^MT-')
rcc2t1[['percent.mt']] <- PercentageFeatureSet(rcc2t1, pattern = '^MT-')
rcc2t2[['percent.mt']] <- PercentageFeatureSet(rcc2t2, pattern = '^MT-')

# Filter out high MT %
#> RCC2N 5089 'cells' -> 2545 cells
#> RCC2T1 6860 'cells' -> 4882 cells
#> RCC2T2 28763 'cells' -> 23900 cells
rcc2n.filt <- subset(rcc2n,
                     subset = #nCount_RNA > 800 &
                       #nFeature_RNA > 500 &
                       percent.mt < 10)
rcc2t1.filt <- subset(rcc2t1,
                      subset = percent.mt < 10)
rcc2t2.filt <- subset(rcc2t2,
                      subset = percent.mt < 10)

###############################################################################|
### Run Doublet Finder
#####
# RCC2N doublet finder; 2545 cells -> 2508 cells
rcc2n.filt$multRate <- 0.016 # from 10X based on cells recovered
rcc2n.filt <- NormalizeData(rcc2n.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc2n@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc2n.filt, PCs = 1:20, sct = F)
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
annotations <- rcc2n.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc2n.filt$multRate)*nrow(rcc2n.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc2n.filt <- doubletFinder(rcc2n.filt,
                            PCs = 1:20,
                            pN = 0.25,
                            pK = pK,
                            nExp = nExp_poi.adj,
                            reuse.pANN = F,
                            sct = F)
colDF <- colnames(rcc2n.filt@meta.data)[grepl('DF.*',colnames(rcc2n.filt@meta.data))]

# Rename doublet identying column
names(rcc2n.filt@meta.data)[names(rcc2n.filt@meta.data) == colDF] <- 'DoubletID'
DimPlot(rcc2n.filt, 
        reduction = 'umap', 
        group.by = 'DoubletID')

# Filter out doublets
rcc2n.filt <- subset(rcc2n.filt, subset = DoubletID == 'Singlet')

# RCC2T1 Doublet Finder; 4882 cells -> 4703 cells
rcc2t1.filt$multRate <- 0.04 # from 10X based on cells recovered
rcc2t1.filt <- NormalizeData(rcc2t1.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc2t1@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc2t1.filt, PCs = 1:20, sct = F)
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
annotations <- rcc2t1.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc2t1.filt$multRate)*nrow(rcc2t1.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc2t1.filt <- doubletFinder(rcc2t1.filt,
                             PCs = 1:20,
                             pN = 0.25,
                             pK = pK,
                             nExp = nExp_poi.adj,
                             reuse.pANN = F,
                             sct = F)
colDF <- colnames(rcc2t1.filt@meta.data)[grepl('DF.*',colnames(rcc2t1.filt@meta.data))]

# Rename doublet identying column
names(rcc2t1.filt@meta.data)[names(rcc2t1.filt@meta.data) == colDF] <- 'DoubletID'
# DimPlot(rcc2t1.filt, 
#         reduction = 'umap', 
#         group.by = 'DoubletID')

# Filter out doublets
rcc2t1.filt <- subset(rcc2t1.filt, subset = DoubletID == 'Singlet')

# RCC2T2 Doublet Finder; 23900 cells -> 22087 cells
rcc2t2.filt$multRate <- 0.08 # from 10X based on cells recovered
rcc2t2.filt <- NormalizeData(rcc2t2.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc2t2@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc2t2.filt, PCs = 1:20, sct = F)
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
annotations <- rcc2t2.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc2t2.filt$multRate)*nrow(rcc2t2.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc2t2.filt <- doubletFinder(rcc2t2.filt,
                             PCs = 1:20,
                             pN = 0.25,
                             pK = pK,
                             nExp = nExp_poi.adj,
                             reuse.pANN = F,
                             sct = F)
colDF <- colnames(rcc2t2.filt@meta.data)[grepl('DF.*',colnames(rcc2t2.filt@meta.data))]

# Rename doublet identying column
names(rcc2t2.filt@meta.data)[names(rcc2t2.filt@meta.data) == colDF] <- 'DoubletID'
# DimPlot(rcc2t2.filt, 
#         reduction = 'umap', 
#         group.by = 'DoubletID')

# Filter out doublets
rcc2t2.filt <- subset(rcc2t2.filt, subset = DoubletID == 'Singlet')
#####

###############################################################################|
### Combine all RCC2 objects
# Merger normal and tissue data
set.seed(555)
rcc2 <- merge(rcc2n.filt, y = c(rcc2t1.filt,rcc2t2.filt),
              add.cell.ids = c("RCC2N", "RCC2T1","RCC2T2")) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc2@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(resolution = 0.8, verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

rcc2$orig.ident <- factor(rcc2$orig.ident,
                          levels = c('RCC2N','RCC2T1','RCC2T2'),
                          labels = c('RCC2N','RCC2T1','RCC2T2'))
table(rcc2$orig.ident)
dimPre <- DimPlot(rcc2, group.by = 'seurat_clusters', split.by = 'orig.ident', label = T)


# Violin plots
violinPre <- VlnPlot(rcc2,
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3,
                     group.by = 'orig.ident')

# Scatter plots
plot1 <- FeatureScatter(rcc2, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident')
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')
plot2

p1o <- DimPlot(object = rcc2, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc2, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

p3o <- DimPlot(object = rcc2,
               reduction = 'umap',
               pt.size = .1,
               group.by = "orig.ident")
# ElbowPlot(rcc2)
qcPlots <- ggarrange(plot1,plot2,p1o,p3o,ncol=2,nrow=2)

# Visualization
orig.stim <- DimPlot(rcc2, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim

condition <- DimPlot(rcc2, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 3,order=T) + ggtitle("By Tumor Normal")
condition
clusters <- DimPlot(rcc2, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 3,order=T) + ggtitle("By UMAP Clusters")
clusters

clusterOriginPlot <- ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

###############################################################################|
### Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc2) <- "RNA" # SCT
rcc2 <- ScaleData(rcc2,verbose = T,features = rownames(rcc2))
es.max = sctype_score(scRNAseqData = rcc2@assays$RNA$scale.data, #SCT
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc2@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc2@meta.data[rcc2@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc2@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc2@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc2@meta.data$cellassign[rcc2@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
# rcc2$ca9 <- rcc2@assays$RNA$scale.data['CA9',]
# rcc2$cellassign <- ifelse(rcc2$cellassign == 'Proximal tubular cell' & rcc2$ca9 > 0,'Proximal Tubular cell + CA9',rcc2$cellassign)
rcc2@meta.data$cellassign <- ifelse(rcc2@meta.data$cellassign == 'Tumor',
                                    paste0(rcc2@meta.data$cellassign,' ',rcc2@meta.data$seurat_clusters),
                                    rcc2@meta.data$cellassign)

# Get gene list of cell types for heatmap
gene.list <- c()
for(i in 1:nrow(kid.mrkrs)){
  gene.list <- c(gene.list, kid.mrkrs[i,]$Symbol)
}
gene.list <- strsplit(gene.list,',')
gene.list <- unlist(gene.list)  

# Plot heatmap
feature_heatmap <- DoHeatmap(rcc2,
                             features = gene.list,#VariableFeatures(rcc2)[1:150],c('CA9','NDUFA4L2','NNMT','VEGFA','HIF1A'),#
                             #cells = 1:500,
                             group.by = 'cellassign',#'seurat_clusters',
                             size = 4,
                             angle = 90) +
  #scale_y_discrete(breaks = c(10,5,5,6,8,8,5,5,6,5,4,4,6,6,10,9,5,8,4,8,5,8,9,4,5,9)) +
  ggtitle('RCC2 Heatmap Tissue')

# Feature plot
cancer.features <- FeaturePlot(rcc2,
                               features = 'CA9',#c('NDUFA4L2','CA9','VEGFA','EGFR','NNMT'),
                               reduction = 'umap',
                               label = T,
                               repel = T,
                               order = T,
                               min.cutoff = 'q10',
                               max.cutoff = 'q90',
                               split.by = 'orig.ident')
cancer.features

custom1 <- DimPlot(rcc2,
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

clusters <- DimPlot(rcc2,
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

custom2 <- DimPlot(rcc2,
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
rcc2.join <- JoinLayers(rcc2)
rcc2.allMrkrs <- FindAllMarkers(rcc2.join,
                                min.pct = 0.25,
                                min.diff.pct = 0.25,
                                verbose = F)
top10 <- rcc2.allMrkrs %>% group_by(cluster) %>% top_n(-10, p_val_adj)
# split dataframe into list if you find that convenient
top10.cids <- split(top10$gene, top10$cluster)

### Cell Percentage Plots
pt <- table(rcc2$cellassign, rcc2$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

colourCount = length(unique(rcc2$cellassign))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cellassign <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('RCC2 Percent Cellassign Composition of Samples')

pt2 <- table(rcc2$seurat_clusters, rcc2$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

colourCount = length(unique(rcc2$seurat_clusters))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cluster <- ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('RCC2 Percent Cluster Composition of Samples')

### ssGSEA workflow
# create column for orig.ident and cellassign
rcc2$orig.cellassign <- paste0(rcc2$orig.ident,'.',rcc2$cellassign)

# Get average gene expression by orig cellassign
aggExpr <- AggregateExpression(rcc2,
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
             main = 'RCC2 Pathways by Aggregated Identity and Cellassign')

### Find gene differences between Normal and Tumor clusters/cell types
# Tumor 1 vs Normal
t1_n <- FindMarkers(rcc2.join,
                    group.by = 'orig.ident',
                    ident.1 = 'RCC2T1',
                    ident.2 = 'RCC2N',
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
t2_n <- FindMarkers(rcc2.join,
                    group.by = 'orig.ident',
                    ident.1 = 'RCC2T2',
                    ident.2 = 'RCC2N',
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

# Tumor 1 vs Normal
t1_t2 <- FindMarkers(rcc2.join,
                     group.by = 'orig.ident',
                     ident.1 = 'RCC2T1',
                     ident.2 = 'RCC2T2',
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

### Investigate tumor cluster differences in T1 and T2
# Create column for origin tissue and cluster
rcc1.join$orig.cluster <- paste0(rcc1$orig.ident,".",rcc1$seurat_clusters)

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