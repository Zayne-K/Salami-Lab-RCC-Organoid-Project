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

#####
### Load in data
# Read in normal tissue
rcc9n.data <- Read10X(data.dir = '/avatar_data4/11559-TP/10x_analysis_11559-TP/Sample_11559-TP-1/filtered_feature_bc_matrix')
rcc9n <- CreateSeuratObject(counts = rcc9n.data, project = 'RCC9N Tissue')

# Read in tumor tissue
rcc9t.data <- Read10X(data.dir = '/avatar_data4/11559-TP/10x_analysis_11559-TP/Sample_11559-TP-2/filtered_feature_bc_matrix/')
rcc9t <- CreateSeuratObject(counts = rcc9t.data, project = 'RCC9T Tissue')

### Perform Mito Filtering
# Get mitochondrial RNA percentage
rcc9n[['percent.mt']] <- PercentageFeatureSet(rcc9n, pattern = '^MT-')
rcc9t[['percent.mt']] <- PercentageFeatureSet(rcc9t, pattern = '^MT-')

# Filter out high MT %
#> RCC9N 8928 'cells' -> 7274 cells
#> RCC9T 5067 'cells' -> 3916 cells
rcc9n.filt <- subset(rcc9n,
                     subset = #nCount_RNA > 800 &
                       #nFeature_RNA > 500 &
                       percent.mt < 10)
rcc9t.filt <- subset(rcc9t,
                     subset = percent.mt < 10)

###############################################################################|
### Run Doublet Finder
# RCC9N doublet finder; 7274 cells -> 7232 cells
rcc9n.filt$multRate <- 0.0061 # from 10X based on cells recovered
rcc9n.filt <- NormalizeData(rcc9n.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc9n@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc9n.filt, PCs = 1:20, sct = F)
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
annotations <- rcc9n.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc9n.filt$multRate)*nrow(rcc9n.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc9n.filt <- doubletFinder(rcc9n.filt,
                            PCs = 1:20,
                            pN = 0.25,
                            pK = pK,
                            nExp = nExp_poi.adj,
                            reuse.pANN = F,
                            sct = F)
colDF <- colnames(rcc9n.filt@meta.data)[grepl('DF.*',colnames(rcc9n.filt@meta.data))]

# Rename doublet identying column
names(rcc9n.filt@meta.data)[names(rcc9n.filt@meta.data) == colDF] <- 'DoubletID'
DimPlot(rcc9n.filt, 
        reduction = 'umap', 
        group.by = 'DoubletID')

# Filter out doublets
rcc9n.filt <- subset(rcc9n.filt, subset = DoubletID == 'Singlet')

# RCC9T Doublet Finder; 3916 cells -> 3903 cells
rcc9t.filt$multRate <- 0.0039 # from 10X based on cells recovered
rcc9t.filt <- NormalizeData(rcc9t.filt, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc9t@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc9t.filt, PCs = 1:20, sct = F)
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
annotations <- rcc9t.filt@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc9t.filt$multRate)*nrow(rcc9t.filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc9t.filt <- doubletFinder(rcc9t.filt,
                            PCs = 1:20,
                            pN = 0.25,
                            pK = pK,
                            nExp = nExp_poi.adj,
                            reuse.pANN = F,
                            sct = F)
colDF <- colnames(rcc9t.filt@meta.data)[grepl('DF.*',colnames(rcc9t.filt@meta.data))]

# Rename doublet identying column
names(rcc9t.filt@meta.data)[names(rcc9t.filt@meta.data) == colDF] <- 'DoubletID'
# DimPlot(rcc9t.filt, 
#         reduction = 'umap', 
#         group.by = 'DoubletID')

# Filter out doublets
rcc9t.filt <- subset(rcc9t.filt, subset = DoubletID == 'Singlet')

###############################################################################|
### Combine all RCC9 objects
# Merger normal and tissue data
set.seed(555)
rcc9 <- merge(rcc9n.filt, y = c(rcc9t.filt),
              add.cell.ids = c("RCC9N", "RCC9T")) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc9@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(resolution = 0.8, verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

rcc9$orig.ident <- factor(rcc9$orig.ident,
                          levels = c('RCC9N Tissue','RCC9T Tissue'),
                          labels = c('RCC9N','RCC9T'))
table(rcc9$orig.ident)
dimPre <- DimPlot(rcc9, group.by = 'seurat_clusters', split.by = 'orig.ident', label = T)


# Violin plots
violinPre <- VlnPlot(rcc9,
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3,
                     group.by = 'orig.ident')

# Scatter plots
plot1 <- FeatureScatter(rcc9, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident')
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')
plot2

p1o <- DimPlot(object = rcc9, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc9, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

p3o <- DimPlot(object = rcc9,
               reduction = 'umap',
               pt.size = .1,
               group.by = "orig.ident")
# ElbowPlot(rcc9)
qcPlots <- ggarrange(plot1,plot2,p1o,p3o,ncol=2,nrow=2)

# Visualization
orig.stim <- DimPlot(rcc9, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim

# Get gene list of cell types for heatmap
gene.list <- c()
for(i in 1:nrow(kid.mrkrs)){
  gene.list <- c(gene.list, kid.mrkrs[i,]$Symbol)
}
gene.list <- strsplit(gene.list,',')
gene.list <- unlist(gene.list)  

# Plot heatmap
feature_heatmap <- DoHeatmap(rcc9,
          features = gene.list,#VariableFeatures(rcc9)[1:150],c('CA9','NDUFA4L2','NNMT','VEGFA','HIF1A'),#
          #cells = 1:500,
          group.by = 'cellassign',#'seurat_clusters',
          size = 4,
          angle = 90) +
  #scale_y_discrete(breaks = c(10,5,5,6,8,8,5,5,6,5,4,4,6,6,10,9,5,8,4,8,5,8,9,4,5,9)) +
  ggtitle('RCC9 Heatmap Tissue')

 # NNMT Feature plot
cancer.features <- FeaturePlot(rcc9,
                               features = 'CA9',#c('NDUFA4L2','CA9','VEGFA','EGFR','NNMT'),
                               reduction = 'umap',
                               label = T,
                               repel = T,
                               order = T,
                               min.cutoff = 'q10',
                               max.cutoff = 'q90',
                               split.by = 'orig.ident')
cancer.features

condition <- DimPlot(rcc9, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 3,order=T) + ggtitle("By Tumor Normal")
condition
clusters <- DimPlot(rcc9, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 3,order=T) + ggtitle("By UMAP Clusters")
clusters

clusterOriginPlot <- ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

###############################################################################|
### Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc9) <- "RNA" # SCT
rcc9 <- ScaleData(rcc9,verbose = T,features = rownames(rcc9))
es.max = sctype_score(scRNAseqData = rcc9@assays$RNA$scale.data, #SCT
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc9@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc9@meta.data[rcc9@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc9@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc9@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc9@meta.data$cellassign[rcc9@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
# rcc9$ca9 <- rcc9@assays$RNA$scale.data['CA9',]
# rcc9$cellassign <- ifelse(rcc9$cellassign == 'Proximal tubular cell' & rcc9$ca9 > 0,'Proximal Tubular cell + CA9',rcc9$cellassign)

custom1 <- DimPlot(rcc9,
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

clusters <- DimPlot(rcc9,
                    reduction = "umap",
                    pt.size = 0.5,
                    split.by='orig.ident',
                    group.by = "seurat_clusters",
                    label = T,
                    repel = T,
                    label.size = 3,
                    order=T) + 
  ggtitle("By UMAP Clusters")

clusterCellassignPlot <- clusters / custom1

custom2 <- DimPlot(rcc9,
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
rcc9.join <- JoinLayers(rcc9)
rcc9.allMrkrs <- FindAllMarkers(rcc9.join,
                                min.pct = 0.25,
                                min.diff.pct = 0.25,
                                verbose = F)
top10 <- rcc9.allMrkrs %>% group_by(cluster) %>% top_n(-10, p_val_adj)
# split dataframe into list if you find that convenient
top10.cids <- split(top10$gene, top10$cluster)

### Cell Percentage Plots
rcc9@meta.data$cellassign <- ifelse(rcc9@meta.data$cellassign == 'Tumor',
                                    paste0(rcc9@meta.data$cellassign,' ',rcc9@meta.data$seurat_clusters),
                                    rcc9@meta.data$cellassign)
pt <- table(rcc9$cellassign, rcc9$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

colourCount = length(unique(rcc9$cellassign))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cellassign <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('RCC9 Percent Cellassign Composition of Samples')

pt2 <- table(rcc9$seurat_clusters, rcc9$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

colourCount = length(unique(rcc9$seurat_clusters))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cluster <- ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('RCC9 Percent Cluster Composition of Samples')

### Get average gene expression by seurat cluster
geneExpCluster <- AggregateExpression(rcc9,
                                      group.by = 'seurat_clusters',
                                      normalization.method = 'LogNormalize',
                                      return.seurat = T,
                                      verbose = F)
#> returns a matrix of logNormalized summed counts by group
geneCluster <- as.data.frame(geneExpCluster@assays$RNA$scale.data)
geneCluster <- as.data.frame(geneCluster)

### Find gene differences between Normal and Tumor clusters/cell types
# Create column for origin tissue and cluster
rcc9.join$orig.cluster <- paste0(rcc9$orig.ident,".",rcc9$seurat_clusters)
# Tumor cluster 0 vs Normal cluster 6
t0_n6 <- FindMarkers(rcc9.join,
                     group.by = 'orig.cluster',
                     ident.1 = 'RCC9T.0',
                     ident.2 = 'RCC9N.6',
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

# Tumor cluster 1 vs Normal cluster 6
t1_n6 <- FindMarkers(rcc9.join,
                     group.by = 'orig.cluster',
                     ident.1 = 'RCC9T.1',
                     ident.2 = 'RCC9N.6',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
# Tumor cluster 2 vs Normal CLuster 6
t2_n6 <- FindMarkers(rcc9.join,
                     group.by = 'orig.cluster',
                     ident.1 = 'RCC9T.2',
                     ident.2 = 'RCC9N.6',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
# Tumor vs Normal cluster 27
t27_n27 <- FindMarkers(rcc9.join,
                       group.by = 'orig.cluster',
                       ident.1 = 'RCC9T.27',
                       ident.2 = 'RCC9N.27',
                       min.pct = 0.25,
                       min.diff.pct = 0.25,
                       verbose = F)
EnhancedVolcano(t27_n27,
                lab = rownames(t27_n27),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Normal Prox-Tub <--Log2FC--> Up in Tumor Prox-Tub',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Proximal Tubule-2 Cells Tumor v Normal',
                pCutoff = 0.05,
                FCcutoff = 1.5)
# Tumor 1 vs Tumor 2
t1_t2 <- FindMarkers(rcc9.join,
                     group.by = 'orig.cluster',
                     ident.1 = 'RCC9T.1',
                     ident.2 = 'RCC9T.2',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
EnhancedVolcano(t1_t2,
                lab = rownames(t1_t2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Tumor Cluster 2 <--Log2FC--> Up in Tumor Cluster 1',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor cells in Cluster 1 vs Cluster 2',
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

# Tumor 1 vs Tumor 0
t1_t0 <- FindMarkers(rcc9.join,
                     group.by = 'orig.cluster',
                     ident.1 = 'RCC9T.1',
                     ident.2 = 'RCC9T.0',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
EnhancedVolcano(t1_t0,
                lab = rownames(t1_t0),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Tumor Cluster 0 <--Log2FC--> Up in Tumor Cluster 1',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor cells in Cluster 1 vs Cluster 0',
                pCutoff = 0.05,
                FCcutoff = 1.5)
t1_t0$genes <- rownames(t1_t0)

t1_t0.mrkrs <- t1_t0 %>% arrange(desc(avg_log2FC))
fold_changes <- t1_t0.mrkrs$avg_log2FC
names(fold_changes) <-t1_t0.mrkrs$genes
t1_t0.gsea <- fgsea(pathways = hm.sym,
                    stats = fold_changes,
                    eps = 0.0,
                    minSize = 15,
                    maxSize = 500)
t1_t0.gsea$comp <- 'Tumor 1 v Tumor 0'

# Tumor 2 vs Tumor 0
t2_t0 <- FindMarkers(rcc9.join,
                     group.by = 'orig.cluster',
                     ident.1 = 'RCC9T.2',
                     ident.2 = 'RCC9T.0',
                     min.pct = 0.25,
                     min.diff.pct = 0.25,
                     verbose = F)
EnhancedVolcano(t2_t0,
                lab = rownames(t2_t0),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Tumor Cluster 0 <--Log2FC--> Up in Tumor Cluster 2',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor cells in Cluster 2 vs Cluster 0',
                pCutoff = 0.05,
                FCcutoff = 1.5)
t2_t0$genes <- rownames(t2_t0)

t2_t0.mrkrs <- t2_t0 %>% arrange(desc(avg_log2FC))
fold_changes <- t2_t0.mrkrs$avg_log2FC
names(fold_changes) <-t2_t0.mrkrs$genes
t2_t0.gsea <- fgsea(pathways = hm.sym,
                    stats = fold_changes,
                    eps = 0.0,
                    minSize = 15,
                    maxSize = 500)
t2_t0.gsea$comp <- 'Tumor 2 v Tumor 0'

# Tumor endothelial vs Normal endothelial
tEndo_nEndo <- FindMarkers(rcc9.join,
                           group.by = 'cellassign',
                           ident.1 = 'RCC9T.0',
                           ident.2 = 'RCC9N.6',
                           min.pct = 0.25,
                           min.diff.pct = 0.25,
                           verbose = F)

###############################################################################|
### Combined FGSEA bubble plot
#####
tops <- rbind(t1_t0.gsea,t2_t0.gsea,t1_t2.gsea)
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
tops$Direction <- factor(ifelse(tops$padj < 0.05 & tops$NES>0,"Enriched in T1/T2",
                                ifelse(tops$padj < 0.05 & tops$NES<0,"Enriched in T0/T2","Not significant")),
                         levels = c("Enriched in T1/T2", "Enriched in T0/T2","Not significant"))
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