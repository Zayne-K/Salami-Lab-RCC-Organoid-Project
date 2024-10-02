rm(list=ls())

.libPaths("/home/kzayne/R/x86_64-pc-linux-gnu-library/R-4.2.3/")

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

# rcc1.aaron <- readRDS("/mnt/DATA2/ccRCC_ssRNAseq/udager_all.rds")

# Load RCC1 datasets
# pathways: .67 = '/mnt/DATA2/ccRCC... ; .80 = /avatar_Data2/ccRCC...
rcc1n.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/sample1Nnew/outs/filtered_feature_bc_matrix")
rcc1n <- CreateSeuratObject(counts = rcc1n.data,project = "RCC1N")
rcc1t1.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/sample1T1new/outs/filtered_feature_bc_matrix")
rcc1t1 <- CreateSeuratObject(counts = rcc1t1.data,project = "RCC1T1")
rcc1t2.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/sample1T2new/outs/filtered_feature_bc_matrix")
rcc1t2 <- CreateSeuratObject(counts = rcc1t2.data,project = "RCC1T2")

# Load RCC2 datasets
rcc2n.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/10x_analysis_5964-SS/Sample_5964-SS-1_ACGACTACCA-TTAGGGTCGT/filtered_feature_bc_matrix")
rcc2n <- CreateSeuratObject(counts = rcc2n.data,project = "RCC2N")
rcc2t1.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/10x_analysis_5964-SS/Sample_5964-SS-2_CCCAGCTTCT-GTTTGGTGTC/filtered_feature_bc_matrix")
rcc2t1 <- CreateSeuratObject(counts = rcc2t1.data,project = "RCC2T1")
rcc2t2.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/10x_analysis_5964-SS/Sample_5964-SS-3_TTGAGAGTCA-CTACCAGGTT/filtered_feature_bc_matrix")
rcc2t2 <- CreateSeuratObject(counts = rcc2t2.data,project = "RCC2T2")

# Load RCC3 datasets
rcc3t1.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/RCC3/SI_32602/filtered_feature_bc_matrix")
rcc3t1 <- CreateSeuratObject(counts = rcc3t1.data,project = "RCC3T1")
rcc3t2.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/RCC3/SI_32603/filtered_feature_bc_matrix")
rcc3t2 <- CreateSeuratObject(counts = rcc3t2.data,project = "RCC3T2")
rcc3n.data <- Read10X(data.dir = "/avatar_data2/ccRCC_ssRNAseq/RCC3/SI_32604/filtered_feature_bc_matrix")
rcc3n <- CreateSeuratObject(counts = rcc3n.data,project = "RCC3N")

rcc3t.org.pre.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC3T_org_nerat_treat/outs/filtered_feature_bc_matrix')
rcc3t.org.pre <- CreateSeuratObject(counts = rcc3t.org.pre.data, project = 'RCC3T Org PreTreat')

rcc3t.org.post.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC3T_org_nerat_treat_no/outs/filtered_feature_bc_matrix')
rcc3t.org.post <- CreateSeuratObject(counts = rcc3t.org.post.data, project = 'RCC3T Org PostTreat')

# Load RCC4 datasets
rcc4t2.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/RCC4/RCC4T2/outs/filtered_feature_bc_matrix')
rcc4t2 <- CreateSeuratObject(counts = rcc4t2.data, project = 'RCC4T2')

# Load RCC5 datasets
rcc5n.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC5N_tissue/outs/filtered_feature_bc_matrix')
rcc5n <- CreateSeuratObject(counts = rcc5n.data, project = 'RCC5N')

rcc5t1.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC5T1_tissue/outs/filtered_feature_bc_matrix')
rcc5t1 <- CreateSeuratObject(counts = rcc5t1.data, project = 'RCC5T1')

rcc5n_org.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC5N_org/outs/filtered_feature_bc_matrix')
rcc5n_org <- CreateSeuratObject(counts = rcc5n_org.data, project = 'RCC5N_org')

rcc5t1t2_org_pre.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC5T1T2_org_nerat_treat_no/outs/filtered_feature_bc_matrix')
rcc5t1t2_org_preTreat <- CreateSeuratObject(counts = rcc5t1t2_org_pre.data, project = 'RCC5_org Pre_Treat')

rcc5t1t2_org_post.data <- Read10X(data.dir = '/avatar_data2/ccRCC_ssRNAseq/10716TP/RCC5T1T2_org_nerat_treat/outs/filtered_feature_bc_matrix')
rcc5t1t2_org_postTreat <- CreateSeuratObject(counts = rcc5t1t2_org_post.data, project = 'RCC5_org Post_Treat')

# Load RCC9 datasets
rcc9n.data <- Read10X(data.dir = '/avatar_data4/11559-TP/10x_analysis_11559-TP/Sample_11559-TP-1/filtered_feature_bc_matrix')
rcc9n <- CreateSeuratObject(counts = rcc9n.data, project = 'RCC9N Tissue')

# Load RCC10 datasets
rcc10n.data <- Read10X(data.dir = '/avatar_data4/11632-TP/10x_analysis_11632-TP/Sample_11632-TP-1/filtered_feature_bc_matrix')
rcc10n <- CreateSeuratObject(counts = rcc10n.data, project = 'RCC10N Tissue')

rcc10t1.data <- Read10X(data.dir = '/avatar_data4/11632-TP/10x_analysis_11632-TP/Sample_11632-TP-2/filtered_feature_bc_matrix/')
rcc10t1 <- CreateSeuratObject(counts = rcc10t1.data, project = 'RCC10T1 Tissue')

rcc10t2.data <- Read10X(data.dir = '/avatar_data4/11632-TP/10x_analysis_11632-TP/Sample_11632-TP-3/filtered_feature_bc_matrix/')
rcc10t2 <- CreateSeuratObject(counts = rcc10t2.data, project = 'RCC10T2 Tissue')

# Combine all seurat objects for RCC1
#####
set.seed(555)
rcc1 <- merge(rcc1n, y = c(rcc1t1,rcc1t2), add.cell.ids = c("RCC1N", "RCC1T1", "RCC1T2")) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc1@var.genes, npcs = 20, verbose = FALSE)

rcc1[["percent.mt"]] <- PercentageFeatureSet(rcc1, pattern = "^MT-")
table(rcc1$orig.ident)

# Violin plots
VlnPlot(rcc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plots
plot1 <- FeatureScatter(rcc1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

p1o <- DimPlot(object = rcc1, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc1, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

set.seed(555)
rcc1 <- rcc1 %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

p1n <- DimPlot(object = rcc1, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2n <- VlnPlot(object = rcc1, features = "harmony_1",group.by = "orig.ident", pt.size = .1)

p1o|p1n
p2o|p2n

set.seed(555)
rcc1 <- rcc1 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

# Visualization
orig.stim <- DimPlot(rcc1, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim

condition <- DimPlot(rcc1, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 5,order=T) + ggtitle("By Tumor Normal")
condition
clusters <- DimPlot(rcc1, reduction = "umap",pt.size = 0.5,split.by = 'orig.ident',group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
clusters

ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

# Identify cluster marker genes
rcc1n.sub <- subset(rcc1,orig.ident == 'RCC1N')
rcc1t1.sub <- subset(rcc1,orig.ident == 'RCC1T1')
rcc1t2.sub <- subset(rcc1,orig.ident == 'RCC1T2')

rcc1n.mrkrs <- FindAllMarkers(rcc1n.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
rcc1n.mrkr.table <- rcc1n.mrkrs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

rcc1t1.mrkrs <- FindAllMarkers(rcc1t1.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
rcc1t1.mrkr.table <- rcc1t1.mrkrs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

rcc1t2.mrkrs <- FindAllMarkers(rcc1t2.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
rcc1t2.mrkr.table <- rcc1t2.mrkrs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# rcc1.mrkrs <-  FindAllMarkers(rcc1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# rcc1.mrkr.table <- rcc1.mrkrs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# Differentiate markers between 9,10,11
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster9.markers <- FindMarkers(rcc1n.sub, ident.1 = 9, ident.2 = c(10, 11), min.pct = 0.25)
head(cluster9.markers, n = 10)
tail(cluster9.markers, n = 10)

cluster10.markers <- FindMarkers(rcc1n.sub, ident.1 = 10, ident.2 = c(9, 11), min.pct = 0.25)
head(cluster9.markers, n = 10)

cluster11.markers <- FindMarkers(rcc1n.sub, ident.1 = 11, ident.2 = c(9, 10), min.pct = 0.25)
head(cluster9.markers, n = 10)

# Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc1) <- "RNA"
rcc1 <- ScaleData(rcc1,verbose = T,features = rownames(rcc1))
es.max = sctype_score(scRNAseqData = rcc1@assays$RNA$scale.data,
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc1@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc1@meta.data[rcc1@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc1@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc1@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc1@meta.data$cellassign[rcc1@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
# rcc1$ca9 <- rcc1@assays$RNA$scale.data['CA9',]
rcc1$cellassign <- ifelse(rcc1$cellassign == 'Tumor' & rcc1$orig.ident == 'RCC1N','Proximal tubule-3',rcc1$cellassign)

# Cluster and Cellassign UMAP
custom1 <- DimPlot(rcc1, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom1

clusters <- DimPlot(rcc1, reduction = "umap",pt.size = 0.5,split.by = 'orig.ident',group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
clusters

custom2 <- DimPlot(rcc1, reduction = "umap", label = TRUE, repel = TRUE,
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom2

# Cluster and Cellassign Umap of each subset of RCC1 individually
DefaultAssay(rcc1n.sub) <- "RNA"
rcc1n.sub <- ScaleData(rcc1n.sub,verbose = T,features = rownames(rcc1n.sub))
es.max = sctype_score(scRNAseqData = rcc1n.sub@assays$RNA$scale.data,
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls.rcc1n = do.call("rbind", lapply(unique(rcc1n.sub@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc1n.sub@meta.data[rcc1n.sub@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc1n.sub@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.rcc1n = cL_resutls.rcc1n %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc1n.sub@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc1n.sub@meta.data$cellassign[rcc1n.sub@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

custom1.rcc1n <- DimPlot(rcc1n.sub, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
# custom1.rcc1n

clusters.rcc1n <- DimPlot(rcc1n.sub, reduction = "umap",pt.size = 0.5,split.by = 'orig.ident',group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
# clusters.rcc1n
clusters.rcc1n / custom1.rcc1n

#####################
DefaultAssay(rcc1t1.sub) <- "RNA"
rcc1t1.sub <- ScaleData(rcc1t1.sub,verbose = T,features = rownames(rcc1t1.sub))
es.max = sctype_score(scRNAseqData = rcc1t1.sub@assays$RNA$scale.data,
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls.rcc1t1 = do.call("rbind", lapply(unique(rcc1t1.sub@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc1t1.sub@meta.data[rcc1t1.sub@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc1t1.sub@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.rcc1t1 = cL_resutls.rcc1t1 %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc1t1.sub@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc1t1.sub@meta.data$cellassign[rcc1t1.sub@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

custom1.rcc1t1 <- DimPlot(rcc1t1.sub, reduction = "umap", label = TRUE, repel = TRUE,
                         split.by = "orig.ident",
                         group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
# custom1.rcc1t1

clusters.rcc1t1 <- DimPlot(rcc1t1.sub, reduction = "umap",pt.size = 0.5,split.by = 'orig.ident',group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
# clusters.rcc1t1
clusters.rcc1t1 / custom1.rcc1t1

#####################
DefaultAssay(rcc1t2.sub) <- "RNA"
rcc1t2.sub <- ScaleData(rcc1t2.sub,verbose = T,features = rownames(rcc1t2.sub))
es.max = sctype_score(scRNAseqData = rcc1t2.sub@assays$RNA$scale.data,
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls.rcc1t2 = do.call("rbind", lapply(unique(rcc1t2.sub@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc1t2.sub@meta.data[rcc1t2.sub@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc1t2.sub@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.rcc1t2 = cL_resutls.rcc1t2 %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc1t2.sub@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc1t2.sub@meta.data$cellassign[rcc1t2.sub@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

custom1.rcc1t2 <- DimPlot(rcc1t2.sub, reduction = "umap", label = TRUE, repel = TRUE,
                          split.by = "orig.ident",
                          group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
# custom1.rcc1t2

clusters.rcc1t2 <- DimPlot(rcc1t2.sub, reduction = "umap",pt.size = 0.5,split.by = 'orig.ident',group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
# clusters.rcc1t2
clusters.rcc1t2 / custom1.rcc1t2

# Normalize data
rcc1.norm <- NormalizeData(rcc1)

# cairo_pdf(paste0(setdir,"RCC1_post_harmony.pdf"),width = 22,height = 20)
# print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
# dev.off()
# 
# cairo_pdf(paste0(setdir,"RCC1 UMAP cell clusters.pdf"),width = 22,height = 11)
# print(custom1)
# dev.off()

FeaturePlot(rcc1,
            repel=F,order = T,
            features = c("NNMT","CA9","KRT18","IGF2BP3"),
            split.by = "orig.ident",pt.size = 0.5)
#####

# Combine all seurat objects for RCC2
#####
set.seed(555)
rcc2 <- merge(rcc2n, y = c(rcc2t1,rcc2t2), add.cell.ids = c("RCC2N", "RCC2T1", "RCC2T2")) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc2@var.genes, npcs = 20, verbose = FALSE)

rcc2[["percent.mt"]] <- PercentageFeatureSet(rcc2, pattern = "^MT-")
table(rcc2$orig.ident)

# Violin plots
VlnPlot(rcc2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plots
plot1 <- FeatureScatter(rcc2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

p1o <- DimPlot(object = rcc2, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc2, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

set.seed(555)
rcc2 <- rcc2 %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

p1n <- DimPlot(object = rcc2, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2n <- VlnPlot(object = rcc2, features = "harmony_1",group.by = "orig.ident", pt.size = .1)

ggarrange(p1o,p1n,p2o,p2n,nrow=2,ncol=2)

set.seed(555)
rcc2 <- rcc2 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

# Visualization
orig.stim <- DimPlot(rcc2, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim

condition <- DimPlot(rcc2, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 5,order=T) + ggtitle("By Tumor Normal")
condition
clusters <- DimPlot(rcc2, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
clusters

ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

# Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc2) <- "RNA"
rcc2 <- ScaleData(rcc2,verbose = T,features = rownames(rcc2))
es.max = sctype_score(scRNAseqData = rcc2@assays$RNA$scale.data,
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
rcc2$ca9 <- rcc2@assays$RNA$scale.data['CA9',]
rcc2$cellassign <- ifelse(rcc2$cellassign == 'Proximal tubular cell' & rcc2$ca9 > 0,'Proximal Tubular cell + CA9',rcc2$cellassign)

custom1 <- DimPlot(rcc2, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom1

clusters <- DimPlot(rcc2, reduction = "umap",pt.size = 0.5,split.by = 'orig.ident', group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
clusters

custom2 <- DimPlot(rcc2, reduction = "umap", label = TRUE, repel = TRUE,
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom2

# cairo_pdf(paste0(setdir,"RCC2_post_harmony.pdf"),width = 22,height = 20)
# print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
# dev.off()
# 
# cairo_pdf(paste0("RCC2 UMAP cell clusters.pdf"),width = 22,height = 11)
# print(custom1)
# dev.off()

rcc2.features <- FeaturePlot(rcc2,
            repel=F,order = T,
            features = c("NNMT","CA9","KRT18","IGF2BP3"),
            split.by = "orig.ident",pt.size = 0.5)
#####

# Combine all seurat objects for RCC3
#####
set.seed(555)
rcc3 <- merge(rcc3n, y = c(rcc3t1,rcc3t2,rcc3t.org.pre,rcc3t.org.post), add.cell.ids = c("RCC3N", "RCC3T1", "RCC3T2","RCC3T Org PreTreat","RCC3T Org PostTreat")) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc3@var.genes, npcs = 20, verbose = FALSE)

rcc3[["percent.mt"]] <- PercentageFeatureSet(rcc3, pattern = "^MT-")
table(rcc3$orig.ident)

# Violin plots
VlnPlot(rcc3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plots
plot1 <- FeatureScatter(rcc3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

p1o <- DimPlot(object = rcc3, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc3, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

# Filter out cells of poor qc
# rcc3
# rcc3 <- subset(rcc3,
#                subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# rcc3

# Run Harmony
set.seed(555)
rcc3 <- rcc3 %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

p1n <- DimPlot(object = rcc3, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2n <- VlnPlot(object = rcc3, features = "harmony_1",group.by = "orig.ident", pt.size = .1)

ggarrange(p1o,p1n,p2o,p2n,nrow=2,ncol=2)

set.seed(555)
rcc3 <- rcc3 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

# Visualization
orig.stim <- DimPlot(rcc3, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim

condition <- DimPlot(rcc3, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 5,order=T) + ggtitle("By Tumor Normal")
condition
clusters <- DimPlot(rcc3, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
clusters

ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

# Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc3) <- "RNA"
rcc3 <- ScaleData(rcc3,verbose = T,features = rownames(rcc3))
es.max = sctype_score(scRNAseqData = rcc3@assays$RNA$scale.data,
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
rcc3$ca9 <- rcc3@assays$RNA$scale.data['CA9',]
rcc3$cellassign <- ifelse(rcc3$cellassign == 'Proximal tubular cell' & rcc3$ca9 > 0,'Proximal Tubular cell + CA9',rcc3$cellassign)

rcc3$orig.ident <- factor(rcc3$orig.ident,
                          levels = c('RCC3N','RCC3T1','RCC3T2','RCC3T Org PreTreat','RCC3T Org PostTreat'))

custom1 <- DimPlot(rcc3, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom1

clusters <- DimPlot(rcc3, reduction = "umap",pt.size = 0.5,split.by='orig.ident',group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")


custom2 <- DimPlot(rcc3, reduction = "umap", label = TRUE, repel = TRUE,
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom2

# cairo_pdf(paste0(setdir,"RCC3_post_harmony.pdf"),width = 22,height = 20)
# print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
# dev.off()
# 
# cairo_pdf(paste0(setdir,"RCC3 UMAP cell clusters.pdf"),width = 22,height = 11)
# print(custom1)
# dev.off()

rcc3.features <- FeaturePlot(rcc3,
            repel=F,order = T,
            features = c("NNMT","CA9",'KRT18','IGF2BP3'),#"NCAM1","EPCAM","IFG2BP3"),
            split.by = "orig.ident",pt.size = 0.5)


# FindMarkers
rcc3.join <- JoinLayers(rcc3, verbose = F)
rcc3.all.mrkrs <- FindAllMarkers(rcc3.join,
                                 min.pct = 0.25,
                                 min.diff.pct = 0.25,
                                 verbose = F)

# subset to organoids
rcc3.org <- subset(rcc3, subset = orig.ident %in% c('RCC3T Org PreTreat','RCC3T Org PostTreat'))
rcc3.org@meta.data$orig.cluster <- paste0(rcc3.org@meta.data$orig.ident,"_",rcc3.org@meta.data$seurat_clusters)
rcc3.org@meta.data$orig.cell <- paste0(rcc3.org@meta.data$orig.ident,"_",rcc3.org@meta.data$cellassign)
rcc3.org.join <- JoinLayers(rcc3.org, verbose = F)

# Find markers between Pre and Post cluster 7
pre_post.c7 <- FindMarkers(rcc3.org.join,
                           group.by = 'orig.cluster',
                           ident.1 = 'RCC3T Org PostTreat_7',
                           ident.2 = 'RCC3T Org PreTreat_7',
                           min.pct = 0.25,
                           min.diff.pct = 0.25,
                           verbose = F)
pre_post.c7$genes <- rownames(pre_post.c7)

EnhancedVolcano(pre_post.c7,
                lab = rownames(pre_post.c7),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Pre-Treated <--Log2FC--> Up in Post-Treated',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 7 cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c7 <- pre_post.c7 %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c7$avg_log2FC
names(fold_changes) <- pre_post.c7$genes
pre_post.c7.gsea <- fgsea(pathways = hm.sym,
                             stats = fold_changes,
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)
pre_post.c7.gsea$comp <- 'Cellassign Cluster 7 Post v Pre'

# Find markers between Pre and Post cluster 1
pre_post.c1.rcc3 <- FindMarkers(rcc3.org.join,
                           group.by = 'orig.cluster',
                           ident.1 = 'RCC3T Org PostTreat_1',
                           ident.2 = 'RCC3T Org PreTreat_1',
                           min.pct = 0.25,
                           min.diff.pct = 0.25,
                           verbose = F)
pre_post.c1.rcc3$genes <- rownames(pre_post.c1.rcc3)

EnhancedVolcano(pre_post.c1.rcc3,
                lab = rownames(pre_post.c1.rcc3),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Pre-Treated <--Log2FC--> Up in Post-Treated',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 1 cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c1.rcc3 <- pre_post.c1.rcc3 %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c1.rcc3$avg_log2FC
names(fold_changes) <- pre_post.c1.rcc3$genes
pre_post.c1.gsea <- fgsea(pathways = hm.sym,
                          stats = fold_changes,
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)
pre_post.c1.gsea$comp <- 'Cellassign Cluster 1 Post v Pre'

# Find markers between Pre and Post cluster 11
pre_post.c11.rcc3 <- FindMarkers(rcc3.org.join,
                                group.by = 'orig.cluster',
                                ident.1 = 'RCC3T Org PostTreat_11',
                                ident.2 = 'RCC3T Org PreTreat_11',
                                min.pct = 0.25,
                                min.diff.pct = 0.25,
                                verbose = F)
pre_post.c11.rcc3$genes <- rownames(pre_post.c11.rcc3)

EnhancedVolcano(pre_post.c11.rcc3,
                lab = rownames(pre_post.c11.rcc3),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Pre-Treated <--Log2FC--> Up in Post-Treated',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 11 cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c11.rcc3 <- pre_post.c11.rcc3 %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c11.rcc3$avg_log2FC
names(fold_changes) <- pre_post.c11.rcc3$genes
pre_post.c11.gsea <- fgsea(pathways = hm.sym,
                          stats = fold_changes,
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)
pre_post.c11.gsea$comp <- 'Cellassign Cluster 11 Post v Pre'

############## Fgseplots combined - bubbleplot
tops <- rbind(pre_post.c7.gsea,pre_post.c1.gsea,pre_post.c11.gsea)
tops <- tops[tops$padj < 0.05,]
tops <- tops[,c(1:7,9)]
tops <- tops[!grepl("HALLMARK",tops$pathway),]
tops$pathway <- gsub("HALLMARK_","",tops$pathway)
tops$pathway <- gsub("_"," ",tops$pathway)
tops <- tops[!grepl("ESTROGEN",tops$pathway),]
tops <- tops[order(tops$padj,decreasing = F),]
tops$pathway <- factor(tops$pathway,levels = unique(c(tops$pathway)))

tops$comp <- factor(tops$comp)

tops$NES <- as.numeric(round(tops$NES,digits=2))
tops$Direction <- factor(ifelse(tops$padj < 0.05 & tops$NES>0,"Enriched in Post",
                                ifelse(tops$padj < 0.05 & tops$NES<0,"Enriched in Pre","Not significant")),
                         levels = c("Enriched in Post", "Enriched in Pre","Not significant"))
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

# Combine all seurat objects for RCC4
#####
set.seed(555)
# rcc4 <- merge(rcc4n, y = c(rcc4t1,rcc4t2), add.cell.ids = c("RCC4N", "RCC4T1", "RCC4T2")) %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE) %>%
#   RunPCA(pc.genes = rcc4@var.genes, npcs = 20, verbose = FALSE)
rcc4 <- NormalizeData(rcc4t2,verbose = F) %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>%
  ScaleData(verbose = F) %>%
  RunPCA(pc.genes = rcc4@var.genes, npcs = 20, verbose = F)

rcc4[["percent.mt"]] <- PercentageFeatureSet(rcc4, pattern = "^MT-")
table(rcc4$orig.ident)

# Violin plots
VlnPlot(rcc4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plots
plot1 <- FeatureScatter(rcc4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

p1o <- DimPlot(object = rcc4, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc4, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

# set.seed(555)
# rcc4 <- rcc4 %>%
#   RunHarmony("orig.ident", plot_convergence = TRUE)
# 
# p1n <- DimPlot(object = rcc4, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
# p2n <- VlnPlot(object = rcc4, features = "harmony_1",group.by = "orig.ident", pt.size = .1)
# 
# ggarrange(p1o,p1n,p2o,p2n,nrow=2,ncol=2)

set.seed(555)
rcc4 <- rcc4 %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>% # harmony
  FindNeighbors(reduction = "pca", dims = 1:20) %>% # harmony
  FindClusters(resolution = 0.5) %>%
  identity()

# Visualization
orig.stim <- DimPlot(rcc4, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim

# condition <- DimPlot(rcc4, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 5,order=T) + ggtitle("By Tumor Normal")
# condition
clusters <- DimPlot(rcc4, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
clusters

ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

# Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc4) <- "RNA"
rcc4 <- ScaleData(rcc4,verbose = T,features = rownames(rcc4))
es.max = sctype_score(scRNAseqData = rcc4@assays$RNA$scale.data,
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc4@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc4@meta.data[rcc4@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc4@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc4@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc4@meta.data$cellassign[rcc4@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
rcc4$ca9 <- rcc4@assays$RNA$scale.data['CA9',]
rcc4$cellassign <- ifelse(rcc4$cellassign == 'Proximal tubular cell' & rcc4$ca9 > 0,'Proximal Tubular cell + CA9',rcc4$cellassign)

custom1 <- DimPlot(rcc4, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom1

clusters <- DimPlot(rcc4, reduction = "umap",pt.size = 0.5,split.by='orig.ident',group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")


custom2 <- DimPlot(rcc4, reduction = "umap", label = TRUE, repel = TRUE,
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom2

# cairo_pdf(paste0(setdir,"RCC4_post_harmony.pdf"),width = 22,height = 20)
# print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
# dev.off()
# 
# cairo_pdf(paste0(setdir,"RCC4 UMAP cell clusters.pdf"),width = 22,height = 11)
# print(custom1)
# dev.off()

rcc4.features <- FeaturePlot(rcc4,
                             repel=F,order = T,
                             features = c("CA9",'KIT','RHCG','LINC01187','FOXI1'),#"NCAM1","EPCAM","IFG2BP3"),
                             split.by = "orig.ident",pt.size = 0.5)

#####

# Combine all RCC9
#####
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
#####

# Combine all RCC10
#####
# Get mitochondrial RNA percentage
rcc10n[['percent.mt']] <- PercentageFeatureSet(rcc10n, pattern = '^MT-')
rcc10t1[['percent.mt']] <- PercentageFeatureSet(rcc10t1, pattern = '^MT-')
rcc10t2[['percent.mt']] <- PercentageFeatureSet(rcc10t2, pattern = '^MT-')
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
rcc10 <- subset(rcc10.pre,
                subset = DoubletID == 'Singlet')

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

#####

################################## RCC5
#####
##### Perform QC, Filtering, and Doublet finder on RCC5 samples
rcc5n[['percent.mt']] <- PercentageFeatureSet(rcc5n, pattern = '^MT-')
rcc5t1[['percent.mt']] <- PercentageFeatureSet(rcc5t1, pattern = '^MT-')
rcc5n_org[['percent.mt']] <- PercentageFeatureSet(rcc5n_org, pattern = '^MT-')
rcc5t1t2_org_preTreat[['percent.mt']] <- PercentageFeatureSet(rcc5t1t2_org_preTreat, pattern = '^MT-')
rcc5t1t2_org_postTreat[['percent.mt']] <- PercentageFeatureSet(rcc5t1t2_org_postTreat, pattern = '^MT-')

# Filter data
rcc5n <- subset(rcc5n,
                subset = #nCount_RNA > 800 &
                  #nFeature_RNA > 500 &
                  percent.mt < 10)
# rcc5t1 <- subset(rcc5t1,
#                 subset = nCount_RNA > 800 &
#                   nFeature_RNA > 500 &
#                   percent.mt < 10)
# rcc5n_org <- subset(rcc5n_org,
#                 subset = nCount_RNA > 800 &
#                   nFeature_RNA > 500 &
#                   percent.mt < 10)
# rcc5t1t2_org_preTreat <- subset(rcc5t1t2_org_preTreat,
#                 subset = nCount_RNA > 800 &
#                   nFeature_RNA > 500 &
#                   percent.mt < 10)
# rcc5t1t2_org_postTreat <- subset(rcc5t1t2_org_postTreat,
#                 subset = nCount_RNA > 800 &
#                   nFeature_RNA > 500 &
#                   percent.mt < 10)

# add multRate for doublet finder
rcc5n$multRate <- 0.004
rcc5n_org$multRate <- 0.023
rcc5t1$multRate <- 0.008
rcc5t1t2_org_preTreat$multRate <- 0.016
rcc5t1t2_org_postTreat$multRate <- 0.004

# Run standard workflow and doublet finder on each sample
#####
samples <- c(rcc5n,rcc5t1,rcc5n_org,rcc5t1t2_org_preTreat,rcc5t1t2_org_postTreat)

rcc5n <- NormalizeData(rcc5n, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc5@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc5n, PCs = 1:20, sct = F)
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
annotations <- rcc5n@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc5n$multRate)*nrow(rcc5n@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc5n <- doubletFinder(rcc5n,
                        PCs = 1:20,
                        pN = 0.25,
                        pK = pK,
                        nExp = nExp_poi.adj,
                        reuse.pANN = F,
                        sct = F)
colDF <- colnames(rcc5n@meta.data)[grepl('DF.*',colnames(rcc5n@meta.data))]
# DimPlot(rcc5n, reduction = 'umap', group.by = colDF)

# Rename doulbet identying column
names(rcc5n@meta.data)[names(rcc5n@meta.data) == colDF] <- 'DoubletID'

# Filter out doublets
rcc5n <- subset(rcc5n, subset = DoubletID == 'Singlet')


rcc5n_org <- NormalizeData(rcc5n_org, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc5@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc5n_org, PCs = 1:20, sct = F)
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
annotations <- rcc5n_org@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc5n_org$multRate)*nrow(rcc5n_org@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc5n_org <- doubletFinder(rcc5n_org,
                        PCs = 1:20,
                        pN = 0.25,
                        pK = pK,
                        nExp = nExp_poi.adj,
                        reuse.pANN = F,
                        sct = F)
colDF <- colnames(rcc5n_org@meta.data)[grepl('DF.*',colnames(rcc5n_org@meta.data))]
# DimPlot(rcc5n_org, reduction = 'umap', group.by = colDF)

# Rename doulbet identying column
names(rcc5n_org@meta.data)[names(rcc5n_org@meta.data) == colDF] <- 'DoubletID'

# Filter out doublets
rcc5n_org <- subset(rcc5n_org, subset = DoubletID == 'Singlet')

rcc5t1 <- NormalizeData(rcc5t1, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc5@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc5t1, PCs = 1:20, sct = F)
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
annotations <- rcc5t1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc5t1$multRate)*nrow(rcc5t1@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc5t1 <- doubletFinder(rcc5t1,
                        PCs = 1:20,
                        pN = 0.25,
                        pK = pK,
                        nExp = nExp_poi.adj,
                        reuse.pANN = F,
                        sct = F)
colDF <- colnames(rcc5t1@meta.data)[grepl('DF.*',colnames(rcc5t1@meta.data))]
# DimPlot(rcc5t1, reduction = 'umap', group.by = colDF)

# Rename doulbet identying column
names(rcc5t1@meta.data)[names(rcc5t1@meta.data) == colDF] <- 'DoubletID'

# Filter out doublets
rcc5t1 <- subset(rcc5t1, subset = DoubletID == 'Singlet')

rcc5t1t2_org_preTreat <- NormalizeData(rcc5t1t2_org_preTreat, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc5@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc5t1t2_org_preTreat, PCs = 1:20, sct = F)
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
annotations <- rcc5t1t2_org_preTreat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc5t1t2_org_preTreat$multRate)*nrow(rcc5t1t2_org_preTreat@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc5t1t2_org_preTreat <- doubletFinder(rcc5t1t2_org_preTreat,
                        PCs = 1:20,
                        pN = 0.25,
                        pK = pK,
                        nExp = nExp_poi.adj,
                        reuse.pANN = F,
                        sct = F)
colDF <- colnames(rcc5t1t2_org_preTreat@meta.data)[grepl('DF.*',colnames(rcc5t1t2_org_preTreat@meta.data))]
# DimPlot(rcc5t1t2_org_preTreat, reduction = 'umap', group.by = colDF)

# Rename doulbet identying column
names(rcc5t1t2_org_preTreat@meta.data)[names(rcc5t1t2_org_preTreat@meta.data) == colDF] <- 'DoubletID'

# Filter out doublets
rcc5t1t2_org_preTreat <- subset(rcc5t1t2_org_preTreat, subset = DoubletID == 'Singlet')

rcc5t1t2_org_postTreat <- NormalizeData(rcc5t1t2_org_postTreat, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc5@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

# find pK based on no ground truth
sweep.res.list_nsclc <- paramSweep(rcc5t1t2_org_postTreat, PCs = 1:20, sct = F)
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
annotations <- rcc5t1t2_org_postTreat@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(unique(rcc5t1t2_org_postTreat$multRate)*nrow(rcc5t1t2_org_postTreat@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))

# run DoubletFinder
rcc5t1t2_org_postTreat <- doubletFinder(rcc5t1t2_org_postTreat,
                        PCs = 1:20,
                        pN = 0.25,
                        pK = pK,
                        nExp = nExp_poi.adj,
                        reuse.pANN = F,
                        sct = F)
colDF <- colnames(rcc5t1t2_org_postTreat@meta.data)[grepl('DF.*',colnames(rcc5t1t2_org_postTreat@meta.data))]
# DimPlot(rcc5t1t2_org_postTreat, reduction = 'umap', group.by = colDF)

# Rename doulbet identying column
names(rcc5t1t2_org_postTreat@meta.data)[names(rcc5t1t2_org_postTreat@meta.data) == colDF] <- 'DoubletID'

# Filter out doublets
rcc5t1t2_org_postTreat <- subset(rcc5t1t2_org_postTreat, subset = DoubletID == 'Singlet')
#####

# Combine all seurat objects for RCC5
#####
set.seed(555)
rcc5 <- merge(rcc5n_org, y = c(rcc5n,rcc5t1,  rcc5t1t2_org_preTreat, rcc5t1t2_org_postTreat),
              add.cell.ids = c("RCC5N_org", "RCC5N", "RCC5T1", "RCC5_org Pre_Treat","RCC5_org Post_Treat")) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = rcc5@var.genes, npcs = 20, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(resolution = 0.8, verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

  
# rcc5 <- merge(rcc5n_org, y = c(rcc5n,rcc5t1,  rcc5t1t2_org_preTreat, rcc5t1t2_org_postTreat),
#               add.cell.ids = c("RCC5N_org", "RCC5N", "RCC5T1", "RCC5_org Pre_Treat","RCC5_org Post_Treat"))
# rcc5 <- SCTransform(rcc5, vars.to.regress = 'percent.mt', verbose = F)
# rcc5 <- RunPCA(rcc5, verbose = FALSE) #  pc.genes = rcc5@var.genes, npcs = 20,
# rcc5 <- FindNeighbors(rcc5, dims = 1:20, verbose = F)
# rcc5 <- FindClusters(rcc5, verbose = F)
# rcc5 <- RunUMAP(rcc5, dims = 1:20, verbose = F)
#   # Seurat::NormalizeData(verbose = FALSE, normalization.method = 'LogNormalize', scale.factor = 10000) %>%

# rcc5 <- VariableFeatures(rcc5,selection.method = "vst")
# rcc5 <- ScaleData(rcc5, verbose = FALSE, vars.to.regress = 'percent.mt')
rcc5$orig.ident <- factor(rcc5$orig.ident,
                          levels = c('RCC5N_org','RCC5N','RCC5T1','RCC5_org Pre_Treat','RCC5_org Post_Treat'))
table(rcc5$orig.ident)
dimPre <- DimPlot(rcc5, group.by = 'seurat_clusters', split.by = 'orig.ident', label = T)


# Violin plots
violingPre <- VlnPlot(rcc5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plots
plot1 <- FeatureScatter(rcc5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

p1o <- DimPlot(object = rcc5, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc5, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

# set.seed(555)
# rcc5 <- rcc5 %>%
#   RunHarmony("orig.ident", plot_convergence = F,assay.use = 'RNA') #SCT

# p1n <- DimPlot(object = rcc5, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
# p2n <- VlnPlot(object = rcc5, features = "harmony_1",group.by = "orig.ident", pt.size = .1)
# 
# qcPlots <- ggarrange(p1o,p1n,p2o,p2n,nrow=2,ncol=2)
ElbowPlot(rcc5)
dimPost <- DimPlot(rcc5, group.by = 'seurat_clusters', split.by = 'orig.ident', label = T, reduction = 'harmony')

# set.seed(555)
# rcc5 <- rcc5 %>% 
#   FindNeighbors(reduction = "pca", dims = 1:20) %>% # harmony
#   FindClusters(resolution = 0.5) %>%
#   RunUMAP(reduction = "pca", dims = 1:20) %>% # harmony
#   identity()
# rcc5$orig.ident <- factor(rcc5$orig.ident,
#                           levels = c('RCC5N_org','RCC5N','RCC5T1','RCC5_org Pre_Treat','RCC5_org Post_Treat'))

# Run doublet finder
#####
suppressMessages(require(DoubletFinder))
rcc5.split <- SplitObject(rcc5,split.by='orig.ident')

for (i in 1:length(rcc5.split)) {
  # print the sample we are on
  print(paste0("Sample ",i))

  # Pre-process seurat object with standard seurat workflow
  rcc.filt.sample <- NormalizeData(rcc5.split[[i]])
  rcc.filt.sample <- FindVariableFeatures(rcc.filt.sample)
  rcc.filt.sample <- ScaleData(rcc.filt.sample)
  rcc.filt.sample <- RunPCA(rcc.filt.sample, nfeatures.print = 10)

  # Find significant PCs
  stdv <- rcc.filt.sample[["pca"]]@stdev
  sum.stdv <- sum(rcc.filt.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                       percent.stdv[2:length(percent.stdv)]) > 0.1),
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc

  # finish pre-processing
  rcc.filt.sample <- RunUMAP(rcc.filt.sample, dims = 1:min.pc)
  rcc.filt.sample <- FindNeighbors(object = rcc.filt.sample, dims = 1:min.pc)
  rcc.filt.sample <- FindClusters(object = rcc.filt.sample, resolution = 0.1)

  # pK identification (no ground-truth)
  sweep.list <- paramSweep(rcc.filt.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)

  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

  ## Homotypic doublet proportion estimate
  annotations <- rcc.filt.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp.poi <- round(optimal.pk * nrow(rcc.filt.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

  # run DoubletFinder
  rcc.filt.sample <- doubletFinder(seu = rcc.filt.sample,
                                   PCs = 1:min.pc,
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
  metadata <- rcc.filt.sample@meta.data
  colnames(metadata)[ncol(metadata)] <- "doublet_finder" #ncol matches number of 'orig.ident' for this group
  rcc.filt.sample@meta.data <- metadata

  # # subset and save
  # rcc.filt.singlets <- subset(rcc.filt.sample, doublet_finder == "Singlet")
  # rcc.filt.split[[i]] <- rcc.filt.singlets
  # remove(rcc.filt.singlets)
  rcc5.split[[i]] <- rcc.filt.sample
}

rcc5.DF.all <- merge(x = rcc5.split[[1]],
                    y = c(rcc5.split[[2]], rcc5.split[[3]], rcc5.split[[4]],
                          rcc5.split[[5]]),
                    project = "RCC5",
                    merge.data = T)
rcc5.DF.all
rcc5.DF.join <- JoinLayers(rcc5.DF.all)
rcc5.DF.join <- RunPCA(rcc5.DF.join)
rcc5.DF.join <- RunUMAP(rcc5.DF.join, dims = 1:10, verbose = F)

DF.name = colnames(rcc5.DF.join@meta.data)[grepl("doublet_finder", colnames(rcc5.DF.join@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(rcc5.DF.join, group.by = "orig.ident") + NoAxes(),
                   DimPlot(rcc5.DF.join, group.by = DF.name) + NoAxes())

VlnPlot(rcc5.DF.all, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

rcc5.filt = rcc5.DF.join[,rcc5.DF.join@meta.data$doublet_finder == "Singlet"]
dim(rcc5.filt)
#####

# Violin plots
VlnPlot(rcc5.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plots
plot1 <- FeatureScatter(rcc5.filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc5.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

violinPost <- VlnPlot(rcc5.filt, features = c('nCount_RNA','nFeature_RNA','percent.mt'),group.by = 'orig.ident',pt.size = 0.1)

p1o <- DimPlot(object = rcc5.filt, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc5.filt, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

# set.seed(555)
# rcc5.filt <- rcc5.filt %>%
#   RunHarmony("orig.ident", plot_convergence = TRUE)
# 
# p1n <- DimPlot(object = rcc5.filt, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
# p2n <- VlnPlot(object = rcc5.filt, features = "harmony_1",group.by = "orig.ident", pt.size = .1)
# 
# qcPlots.postFilter <- ggarrange(p1o,p1n,p2o,p2n,nrow=2,ncol=2)
# ElbowPlot(rcc5.filt)

dimDF <- DimPlot(rcc5.filt, group.by = 'seurat_clusters', split.by = 'orig.ident', label = T, reduction = 'pca')

# Visualization
orig.stim <- DimPlot(rcc5.filt, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By seurat clusters")
orig.stim

# NNMT Feature plot
cancer.features <- FeaturePlot(rcc5.filt,
                    features = c('NDUFA4L2','CA9','VEGFA','EGFR'),
                    reduction = 'umap',
                    label = T,
                    repel = T,
                    order = T,
                    split.by = 'orig.ident')
cancer.features

condition <- DimPlot(rcc5.filt, reduction = "umap",pt.size = 0.5,group.by = "orig.ident",label = T,repel = T,label.size = 3,order=T) + ggtitle("By Tumor Normal")
condition
clusters <- DimPlot(rcc5.filt, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 3,order=T) + ggtitle("By UMAP Clusters")
clusters

ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

# Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc5) <- "RNA" # SCT
rcc5 <- ScaleData(rcc5,verbose = T,features = rownames(rcc5))
es.max = sctype_score(scRNAseqData = rcc5@assays$RNA$scale.data, #SCT
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc5@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc5@meta.data[rcc5@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc5@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc5@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc5@meta.data$cellassign[rcc5@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
# rcc5$ca9 <- rcc5@assays$RNA$scale.data['CA9',]
# rcc5$cellassign <- ifelse(rcc5$cellassign == 'Proximal tubular cell' & rcc5$ca9 > 0,'Proximal Tubular cell + CA9',rcc5$cellassign)

custom1 <- DimPlot(rcc5, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',pt.size = 0.5,label.size = 3,order = T) +
  ggtitle("Cellassign")
custom1

clusters <- DimPlot(rcc5,
                    reduction = "umap",
                    pt.size = 0.5,
                    split.by='orig.ident',
                    group.by = "seurat_clusters",
                    label = T,
                    repel = T,
                    label.size = 3,
                    order=T) + 
  ggtitle("By UMAP Clusters")

# rcc5 <- PrepSCTFindMarkers(rcc5,verbose = F)
rcc5.join <- JoinLayers(rcc5)
rcc5.allMrkrs <- FindAllMarkers(rcc5.join,
                                min.pct = 0.25,
                                min.diff.pct = 0.25,
                                verbose = F)
top10 <- rcc5.allMrkrs %>% group_by(cluster) %>% top_n(-10, p_val_adj)
# split dataframe into list if you find that convenient
top10.cids <- split(top10$gene, top10$cluster)
cl.0413 <- top10[top10$cluster %in% c('0','4','13'),]

rcc5.tumor <- subset(rcc5.join,subset = seurat_clusters %in% c('0','4','13'))
FeaturePlot(rcc5.tumor,
            features = c('CA9','NDUFA4L2'),
            split.by = 'seurat_clusters',
            order = T,
            label = T,
            repel = T)
# consMarkers <- FindConservedMarkers(rcc5.tumor,ident.1 = 0)
rcc5.tissue.c2.mrkrs <- FindMarkers(rcc5,
                                    group.by = 'ident_cluster',
                                    ident.1 = 'RCC5T1_2',
                                    ident.2 = 'RCC5N_2',
                                    min.pct=0.25,
                                    min.diff.pct=0.25,
                                    verbose=F)
rcc5.tissue.c2.mrkrs$genes <- rownames(rcc5.tissue.c2.mrkrs)


# Cell Percentage Plots
rcc5@meta.data$cellassign <- ifelse(rcc5@meta.data$cellassign == 'Tumor',
                                    paste0(rcc5@meta.data$cellassign,' ',rcc5@meta.data$seurat_clusters),
                                    rcc5@meta.data$cellassign)
pt <- table(rcc5$cellassign, rcc5$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

colourCount = length(unique(rcc5$cellassign))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cellassign <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('RCC5 Percent Cellassign Composition of Samples')

pt2 <- table(rcc5$seurat_clusters, rcc5$orig.ident)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

colourCount = length(unique(rcc5$seurat_clusters))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


pct.cluster <- ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('RCC5 Percent Cluster Composition of Samples')

# custom2 <- DimPlot(rcc5.orig, reduction = "umap", label = TRUE, repel = TRUE,
#                    split.by = 'orig.ident',
#                    group.by = 'seurat_clusters',pt.size = 0.5,label.size = 3,order = T) +
#   ggtitle("By UMAP Clusters")
# custom2

# cairo_pdf(paste0(setdir,"RCC5_post_harmony.pdf"),width = 22,height = 20)
# print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
# dev.off()
# 
# cairo_pdf(paste0(setdir,"RCC5 UMAP cell clusters.pdf"),width = 22,height = 11)
# print(custom1)
# dev.off()

rcc5.features <- FeaturePlot(rcc5,
                             repel=F,order = T,
                             features = c("EGFR",'CA9'),#"NCAM1","EPCAM","IFG2BP3"),
                             split.by = "orig.ident",pt.size = 0.5)

# Find markers
rcc5.join <- JoinLayers(rcc5)
rcc5.mrkrs <- FindAllMarkers(rcc5.join,
                             min.pct = 0.25,
                             min.diff.pct = 0.25,
                             logfc.threshold = 0.75,
                             return.thresh = 0.05,
                             verbose = F)

### Get average gene expression by seurat cluster
geneExpCluster <- AggregateExpression(rcc5,
                                    group.by = 'seurat_clusters',
                                    normalization.method = 'LogNormalize',
                                    return.seurat = T,
                                    verbose = F)
#> returns a matrix of logNormalized summed counts by group
geneCluster <- as.data.frame(geneExpCluster@assays$RNA$scale.data)
geneCluster <- as.data.frame(geneCluster)
names(geneCluster) <- c('Cluster0','Cluster1','Cluster2','Cluster3','Cluster4','Cluster5','Cluster6','Cluster7','Cluster8','Cluster9',
                        'Cluster10','Cluster11','Cluster12','Cluster13','Cluster14','Cluster15')

# Run gsea on each cluster
#> Initialize where the
# cluster0 <- as.data.frame(geneCluster[,'Cluster0'])
# rownames(cluster0) <- rownames(geneCluster)
# colnames(cluster0) <- 'expression'
# cluster0$genes <- rownames(cluster0)
# ranks <- cluster0[order(cluster0$expression),]
# ranks <- as.matrix(ranks)
# names(ranks) <- c('expression','genes')
cluster <- rcc5.mrkrs[rcc5.mrkrs$cluster == '15',]
cluster <- cluster %>% arrange(desc(avg_log2FC))
ranks <- cluster$avg_log2FC
names(ranks) <- cluster$gene
pathRes15 <- fgsea(hm.sym,
                 ranks,
                 minSize = 15,
                 maxSize = 500)
pathRes15$comp <- 'Cluster 15'

### Pathway bubble plot - note: need to run for all clusters
tops <- rbind(pathRes0,pathRes1,pathRes2,pathRes3,pathRes4,pathRes5,pathRes6,
              pathRes7,pathRes8,pathRes9,pathRes10,pathRes11,pathRes12,pathRes13,
              pathRes14,pathRes15)
tops <- tops[tops$padj < 0.05,]
tops <- tops[,c(1:7,9)]
tops <- tops[!grepl("HALLMARK",tops$pathway),]
tops$pathway <- gsub("HALLMARK_","",tops$pathway)
tops$pathway <- gsub("_"," ",tops$pathway)
# tops <- tops[!grepl("ESTROGEN",tops$pathway),]
tops <- tops[order(tops$padj,decreasing = F),]
tops$pathway <- factor(tops$pathway,levels = unique(c(tops$pathway)))

tops$comp <- factor(tops$comp,
                    levels = c('Cluster 2','Cluster 3','Cluster 5','Cluster 7',
                               'Cluster 8','Cluster 10','Cluster 11','Cluster 13',
                               'Cluster 14'))

tops$NES <- as.numeric(round(tops$NES,digits=2))
tops$Direction <- factor(ifelse(tops$padj < 0.05 & tops$NES>0,"Up",
                                ifelse(tops$padj < 0.05 & tops$NES<0,"Down","Not significant")),
                         levels = c("Up", "Down","Not significant"))
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
  
  ggtitle("Signigicant Hallmark Pathways for each cluster") +
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


# library(multtest)
# library(metap)
# rcc5.filt.conserved.markers <- data.frame()
# for(cluster in unique(rcc5.filt@meta.data$seurat_clusters)){
#   mrkr.table <- FindConservedMarkers(rcc5.filt.join,
#                                      ident.1 = cluster,
#                                      grouping.var = 'orig.ident',
#                                      verbose = F)
#   rcc5.filt.conserved.markers <- rbind(rcc5.filt.conserved.markers,mrkr.table)
# }

##### Investigation of Mito and Doublet Filtered Clusters 0,4,and 13
rcc5.0413 <- subset(rcc5, subset = seurat_clusters %in% c(0,4,13))
DimPlot(rcc5.0413,
        reduction = 'umap',
        group.by = 'cellassign',
        split.by = 'orig.ident',
        label = T)

FeaturePlot(rcc5.0413,
            features = c('NDUFA4L2','CA9','VEGFA','EGFR'),
            reduction = 'umap',
            label = T,
            repel = T,
            order = T,
            split.by = 'orig.ident')

rcc5.0413.join <- JoinLayers(rcc5.0413)
rcc5.0413.mrkrs <- FindAllMarkers(rcc5.0413.join,
                                  # min.pct = 0.25,
                                  # min.diff.pct = 0.25,
                                  only.pos = T,
                                  verbose = F)
rcc5.0413.top.mrkrs <- function(x){
  rcc5.0413.mrkrs[rcc5.0413.mrkrs$cluster == x, ] %>%
    head(n=10)
}
top10_markers <- rcc5.0413.top.mrkrs(c(0,4,13))

### DEG and GSEA of the Pre and Post treated organoids
rcc5.filt.org <-subset(rcc5.0413, subset = orig.ident %in% c('RCC5_org Post_Treat','RCC5_org Pre_Treat')) #rcc5.filt
rcc5.filt.org@meta.data$ident_cluster <- paste0(rcc5.filt.org@meta.data$orig.ident,"_",rcc5.filt.org@meta.data$seurat_clusters)
rcc5.filt.org@meta.data$ident_cell <- paste0(rcc5.filt.org@meta.data$orig.ident,'_',rcc5.filt.org@meta.data$cellassign)
rcc5.filt.org.join <- JoinLayers(rcc5.filt.org)

# Tumor of Post vs Pre-treated
pre_post.tumor.mrkrs <- FindMarkers(rcc5.filt.org.join,
                                  group.by = 'ident_cell',
                                  ident.1 = 'RCC5_org Post_Treat_Tumor',
                                  ident.2 = 'RCC5_org Pre_Treat_Tumor',
                                  min.pct = 0.25,
                                  min.diff.pct = 0.25,
                                  verbose = F)
pre_post.tumor.mrkrs$genes <- rownames(pre_post.tumor.mrkrs)

EnhancedVolcano(pre_post.tumor.mrkrs,
                lab = rownames(pre_post.tumor.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Pre-Treated <--Log2FC--> Up in Post-Treated',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Tumor cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.tumor.mrkrs <- pre_post.tumor.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.tumor.mrkrs$avg_log2FC
names(fold_changes) <- pre_post.tumor.mrkrs$genes
pre_post.tumor.gsea <- fgsea(pathways = hm.sym,
                             stats = fold_changes,
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)
pre_post.tumor.gsea$comp <- 'Cellassign Tumor Post v Pre'

# Pre Post cluster 0
pre_post.c0.mrkrs <- FindMarkers(rcc5.filt.org.join,
                                    group.by = 'ident_cluster',
                                    ident.1 = 'RCC5_org Post_Treat_0',
                                    ident.2 = 'RCC5_org Pre_Treat_0',
                                    min.pct = 0.25,
                                    min.diff.pct = 0.25,
                                    verbose = F)
pre_post.c0.mrkrs$genes <- rownames(pre_post.c0.mrkrs)

EnhancedVolcano(pre_post.c0.mrkrs,
                lab = rownames(pre_post.c0.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Pre-Treated <--Log2FC--> Up in Post-Treated',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 0 tumor cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c0.mrkrs <- pre_post.c0.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c0.mrkrs$avg_log2FC
names(fold_changes) <- pre_post.c0.mrkrs$genes
pre_post.c0.gsea <- fgsea(pathways = hm.sym,
                             stats = fold_changes,
                             eps = 0.0,
                             minSize = 15,
                             maxSize = 500)
pre_post.c0.gsea$comp <- 'Cluster 0 Post v Pre'

# Pre Post cluster 13
pre_post.c13.mrkrs <- FindMarkers(rcc5.filt.org.join,
                                 group.by = 'ident_cluster',
                                 ident.1 = 'RCC5_org Post_Treat_13',
                                 ident.2 = 'RCC5_org Pre_Treat_13',
                                 min.pct = 0.25,
                                 min.diff.pct = 0.25,
                                 verbose = F)
pre_post.c13.mrkrs$genes <- rownames(pre_post.c13.mrkrs)

EnhancedVolcano(pre_post.c13.mrkrs,
                lab = rownames(pre_post.c13.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Pre-Treated <--Log2FC--> Up in Post-Treated',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 13 tumor cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c13.mrkrs <- pre_post.c13.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c13.mrkrs$avg_log2FC
names(fold_changes) <- pre_post.c13.mrkrs$genes
pre_post.c13.gsea <- fgsea(pathways = hm.sym,
                          stats = fold_changes,
                          eps = 0.0,
                          minSize = 15,
                          maxSize = 500)
pre_post.c13.gsea$comp <- 'Cluster 13 Post v Pre'

# Pre Post cluster 4
pre_post.c4.mrkrs <- FindMarkers(rcc5.filt.org.join,
                                  group.by = 'ident_cluster',
                                  ident.1 = 'RCC5_org Post_Treat_4',
                                  ident.2 = 'RCC5_org Pre_Treat_4',
                                  min.pct = 0.25,
                                  min.diff.pct = 0.25,
                                  verbose = F)
pre_post.c4.mrkrs$genes <- rownames(pre_post.c4.mrkrs)

EnhancedVolcano(pre_post.c4.mrkrs,
                lab = rownames(pre_post.c4.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Pre-Treated <--Log2FC--> Up in Post-Treated',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 4 tumor cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c4.mrkrs <- pre_post.c4.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c4.mrkrs$avg_log2FC
names(fold_changes) <- pre_post.c4.mrkrs$genes
pre_post.c4.gsea <- fgsea(pathways = hm.sym,
                           stats = fold_changes,
                           eps = 0.0,
                           minSize = 15,
                           maxSize = 500)
pre_post.c4.gsea$comp <- 'Cluster 4 Post v Pre'

# Cluster 13 vs 0
pre_post.c130.mrkrs <- FindMarkers(rcc5.filt.org.join,
                                  group.by = 'seurat_clusters',
                                  ident.1 = 13,
                                  ident.2 = 0,
                                  min.pct = 0.25,
                                  min.diff.pct = 0.25,
                                  verbose = F)
pre_post.c130.mrkrs$genes <- rownames(pre_post.c130.mrkrs)

EnhancedVolcano(pre_post.c130.mrkrs,
                lab = rownames(pre_post.c130.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Cluster 0 <--Log2FC--> Up in Cluster 13',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 13 vs Cluster 0 tumor cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c130.mrkrs <- pre_post.c130.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c130.mrkrs$avg_log2FC
names(fold_changes) <- pre_post.c130.mrkrs$genes
pre_post.c130.gsea <- fgsea(pathways = hm.sym,
                           stats = fold_changes,
                           eps = 0.0,
                           minSize = 15,
                           maxSize = 500)
pre_post.c130.gsea$comp <- 'Cluster 13 v Cluster 0'

# Cluster 13 vs 0 Pre only
rcc5.filt.org.pre <- subset(rcc5.filt.org.join,subset = orig.ident == 'RCC5_org Pre_Treat')
pre.c130.mrkrs <- FindMarkers(rcc5.filt.org.pre,
                                   group.by = 'seurat_clusters',
                                   ident.1 = 13,
                                   ident.2 = 0,
                                   min.pct = 0.25,
                                   min.diff.pct = 0.25,
                                   verbose = F)
pre.c130.mrkrs$genes <- rownames(pre.c130.mrkrs)

EnhancedVolcano(pre.c130.mrkrs,
                lab = rownames(pre.c130.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Cluster 0 <--Log2FC--> Up in Cluster 13',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 13 vs Cluster 0 tumor cells in Pre-treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre.c130.mrkrs <- pre.c130.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- pre.c130.mrkrs$avg_log2FC
names(fold_changes) <- pre.c130.mrkrs$genes
pre.c130.gsea <- fgsea(pathways = hm.sym,
                            stats = fold_changes,
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)
pre.c130.gsea$comp <- 'Cluster 13 v Cluster 0 Pre-Treat'

# Cluster 13 vs 0 Post-treated
rcc5.filt.org.post <- subset(rcc5.filt.org.join, subset = orig.ident == 'RCC5_org Post_Treat')
post.c130.mrkrs <- FindMarkers(rcc5.filt.org.post,
                                   group.by = 'seurat_clusters',
                                   ident.1 = 13,
                                   ident.2 = 0,
                                   min.pct = 0.25,
                                   min.diff.pct = 0.25,
                                   verbose = F)
post.c130.mrkrs$genes <- rownames(post.c130.mrkrs)

EnhancedVolcano(post.c130.mrkrs,
                lab = rownames(post.c130.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Cluster 0 <--Log2FC--> Up in Cluster 13',
                ylab = 'Adjuste P-value',
                title = 'DEGs in Cluster 13 vs Cluster 0 tumor cells in Post-treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

post.c130.mrkrs <- post.c130.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- post.c130.mrkrs$avg_log2FC
names(fold_changes) <- post.c130.mrkrs$genes
post.c130.gsea <- fgsea(pathways = hm.sym,
                            stats = fold_changes,
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)
post.c130.gsea$comp <- 'Cluster 13 v Cluster 0 Post-Treat'

# Cluster 13 vs 4
pre_post.c134.mrkrs <- FindMarkers(rcc5.filt.org.join,
                                   group.by = 'seurat_clusters',
                                   ident.1 = 13,
                                   ident.2 = 4,
                                   min.pct = 0.25,
                                   min.diff.pct = 0.25,
                                   verbose = F)
pre_post.c134.mrkrs$genes <- rownames(pre_post.c134.mrkrs)

EnhancedVolcano(pre_post.c134.mrkrs,
                lab = rownames(pre_post.c134.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Cluster 4 <--Log2FC--> Up in Cluster 13',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 13 vs Cluster 4 tumor cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c134.mrkrs <- pre_post.c134.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c134.mrkrs$avg_log2FC
names(fold_changes) <- pre_post.c134.mrkrs$genes
pre_post.c134.gsea <- fgsea(pathways = hm.sym,
                            stats = fold_changes,
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)
pre_post.c134.gsea$comp <- 'Cluster 13 v Cluster 4'

# Cluster 0 vs 4
pre_post.c04.mrkrs <- FindMarkers(rcc5.filt.org.join,
                                   group.by = 'seurat_clusters',
                                   ident.1 = 0,
                                   ident.2 = 4,
                                   min.pct = 0.25,
                                   min.diff.pct = 0.25,
                                   verbose = F)
pre_post.c04.mrkrs$genes <- rownames(pre_post.c04.mrkrs)

EnhancedVolcano(pre_post.c04.mrkrs,
                lab = rownames(pre_post.c04.mrkrs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = 'Up in Cluster 4 <--Log2FC--> Up in Cluster 0',
                ylab = 'Adjusted P-value',
                title = 'DEGs in Cluster 0 vs Cluster 4 tumor cells in Pre and Post treated Organoids',
                pCutoff = 0.05,
                FCcutoff = 1.5)

pre_post.c04.mrkrs <- pre_post.c04.mrkrs %>% arrange(desc(avg_log2FC))
fold_changes <- pre_post.c04.mrkrs$avg_log2FC
names(fold_changes) <- pre_post.c04.mrkrs$genes
pre_post.c04.gsea <- fgsea(pathways = hm.sym,
                            stats = fold_changes,
                            eps = 0.0,
                            minSize = 15,
                            maxSize = 500)
pre_post.c04.gsea$comp <- 'Cluster 0 v Cluster 4'

############## Fgseplots combined - bubbleplot
# tops <- rbind(pre_post.tumor.gsea,pre_post.c0.gsea,pre_post.c13.gsea,pre_post.c130.gsea,pre.c130.gsea,post.c130.gsea)
tops <- rbind(pre_post.c130.gsea,pre_post.c134.gsea,pre_post.c04.gsea)
tops <- tops[tops$padj < 0.05,]
tops <- tops[,c(1:7,9)]
tops <- tops[!grepl("HALLMARK",tops$pathway),]
tops$pathway <- gsub("HALLMARK_","",tops$pathway)
tops$pathway <- gsub("_"," ",tops$pathway)
tops <- tops[!grepl("ESTROGEN",tops$pathway),]
tops <- tops[order(tops$padj,decreasing = F),]
tops$pathway <- factor(tops$pathway,levels = unique(c(tops$pathway)))

tops$comp <- factor(tops$comp)

tops$NES <- as.numeric(round(tops$NES,digits=2))
tops$Direction <- factor(ifelse(tops$padj < 0.05 & tops$NES>0,"Enriched in 13/0",
                                ifelse(tops$padj < 0.05 & tops$NES<0,"Enriched in 0/4","Not significant")),
                         levels = c("Enriched in Enriched in 13/0", "Enriched in 0/4","Not significant"))
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

# Feature Plot of CA9 for clusters 0, 5, 13 only
rcc5.filt.c1305 <- subset(rcc5.filt.org.join, subset = seurat_clusters %in% c(0,5,13))
clusters.c1305 <- DimPlot(rcc5.filt.c1305,
                          group.by = 'seurat_clusters',
                          split.by = 'orig.ident')
ca9_c130 <- FeaturePlot(rcc5.filt.c1305,
                        features = c('CA9','EGFR','VEGFA'),#CXCL14
                        split.by = 'orig.ident')
angio <- FeaturePlot(rcc5.filt.c1305,
                        features = c('LAP3','IFI44','IFITM2','MX1'),
                        split.by = 'orig.ident')
angio.all <- FeaturePlot(rcc5.filt,
                         features = c('LAP3','IFI44','IFITM2','MX1'),
                         split.by = 'orig.ident')

ggarrange(ncol = 2,nrow = 1,clusters.c1305,ca9_c130)

c1305.ca9 <- data.matrix(rcc5.filt.c1305@assays$RNA$counts[,rownames(rcc5.filt.c1305@meta.data[rcc5.filt.c1305@meta.data$seurat_clusters %in% c(0,5,13)&rcc5.filt.c1305@meta.data$orig.ident == 'RCC5T1T2_org Nerat Treated',])])
c1305.ca9.dge <- DGEList(c1305.ca9,remove.zeros = T)
c1305.ca9.dge <- calcNormFactors(c1305.ca9.dge)
c1305.ca9 <- cpm(c1305.ca9.dge,log = T)
c1305.ca9.t <- as.data.frame(t(c1305.ca9[rownames(c1305.ca9) %in% c('EGFR','CA9'),]))
all(rownames(c1305.ca9.t)==rownames(rcc5.filt.c1305@meta.data[rcc5.filt.c1305@meta.data$seurat_clusters %in% c(0,5,13)&rcc5.filt.c1305@meta.data$orig.ident == 'RCC5T1T2_org Nerat Treated',]))
c1305.ca9.t$cluster <- rcc5.filt.c1305@meta.data[rownames(rcc5.filt.c1305@meta.data)%in%rownames(c1305.ca9.t),]$seurat_clusters
boxplot(c1305.ca9.t$EGFR~c1305.ca9.t$cluster)
hist(c1305.ca9.t$CA9)
#####

# Combine all seurat objects for RCC1 RCC2 RCC3 RCC4T2
#####
set.seed(555)
rcc.all <- merge(rcc1n, y = c(rcc1t1,rcc1t2,rcc2n,rcc2t1,rcc2t2,rcc3n,rcc3t1,rcc3t2,rcc4t2),
                 add.cell.ids = c("RCC1N","RCC1T1","RCC1T2","RCC2N","RCC2T1","RCC2T2","RCC3N","RCC3T1","RCC3T2",'RCC4T2')) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = rcc.all@var.genes, npcs = 20, verbose = FALSE)

rcc.all[["percent.mt"]] <- PercentageFeatureSet(rcc.all, pattern = "^MT-")
rcc.all[["percent.rb"]] <- PercentageFeatureSet(rcc.all, pattern = "^RP[SL]")
table(rcc.all$orig.ident)

rm(rcc1n.data,rcc1n,rcc1t1.data,rcc1t1,rcc1t2.data,rcc1t2,
   rcc2n.data,rcc2n,rcc2t1.data,rcc2t1,rcc2t2.data,rcc2t2,
   rcc3n.data,rcc3n,rcc3t1.data,rcc3t1,rcc3t2.data,rcc3t2,
   rcc1,rcc2,rcc3)
gc()

# Violin plots
ggviolin(rcc.all@meta.data,x="orig.ident",y="nFeature_RNA")
ggviolin(rcc.all@meta.data,x="orig.ident",y="nCount_RNA")
ggviolin(rcc.all@meta.data,x="orig.ident",y="percent.mt")
ggviolin(rcc.all@meta.data,x="orig.ident",y="percent.rb")

# Scatter plots
plot1 <- FeatureScatter(rcc.all, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
plot1

# Scatter plots
plot2 <- FeatureScatter(rcc.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
plot2

# QC filters
selected_c <- WhichCells(rcc.all, expression = nFeature_RNA > 200)
selected_f <- rownames(rcc.all)[Matrix::rowSums(rcc.all) > 3]

rcc.filt <- subset(rcc.all, features = selected_f, cells = selected_c)
dim(rcc.all)
dim(rcc.filt)

# rcc.filt.join <- JoinLayers()

# Run doublet finder
suppressMessages(require(DoubletFinder))
rcc.filt.split <- SplitObject(rcc.filt,split.by='orig.ident')

for (i in 1:length(rcc.filt.split)) {
  # print the sample we are on
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  rcc.filt.sample <- NormalizeData(rcc.filt.split[[i]])
  rcc.filt.sample <- FindVariableFeatures(rcc.filt.sample)
  rcc.filt.sample <- ScaleData(rcc.filt.sample)
  rcc.filt.sample <- RunPCA(rcc.filt.sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- rcc.filt.sample[["pca"]]@stdev
  sum.stdv <- sum(rcc.filt.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  rcc.filt.sample <- RunUMAP(rcc.filt.sample, dims = 1:min.pc)
  rcc.filt.sample <- FindNeighbors(object = rcc.filt.sample, dims = 1:min.pc)              
  rcc.filt.sample <- FindClusters(object = rcc.filt.sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(rcc.filt.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- rcc.filt.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(rcc.filt.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  rcc.filt.sample <- doubletFinder(seu = rcc.filt.sample, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
  metadata <- rcc.filt.sample@meta.data
  colnames(metadata)[ncol(metadata)] <- "doublet_finder" #ncol matches number of 'orig.ident' for this group
  rcc.filt.sample@meta.data <- metadata 
  
  # # subset and save
  # rcc.filt.singlets <- subset(rcc.filt.sample, doublet_finder == "Singlet")
  # rcc.filt.split[[i]] <- rcc.filt.singlets
  # remove(rcc.filt.singlets)
  rcc.filt.split[[i]] <- rcc.filt.sample
}

# converge rcc.filt.split
rcc.DF.all <- merge(x = rcc.filt.split[[1]],
                        y = c(rcc.filt.split[[2]], rcc.filt.split[[3]], rcc.filt.split[[4]],
                              rcc.filt.split[[5]], rcc.filt.split[[6]], rcc.filt.split[[7]],
                              rcc.filt.split[[8]],rcc.filt.split[[9]], rcc.filt.split[[10]]),
                        project = "rcc.filt scRNAseq",
                    merge.data = T)
rcc.DF.all
rcc.DF.join <- JoinLayers(rcc.DF.all)
rcc.DF.join <- RunPCA(rcc.DF.join)
rcc.DF.join <- RunUMAP(rcc.DF.join, dims = 1:10, verbose = F)

# rcc.filt = FindVariableFeatures(rcc.filt, verbose = F)
# rcc.filt = ScaleData(rcc.filt, vars.to.regress = c("nFeature_RNA", "percent.mt"),
#                      verbose = F)
# rcc.filt = RunPCA(rcc.filt, verbose = F, npcs = 20)
# rcc.filt = RunUMAP(rcc.filt, dims = 1:10, verbose = F)
# 
# # define the expected number of doublet cells.
# sweep.res.list_kidney <- paramSweep(rcc.filt, PCs = 1:10, sct = FALSE)
# sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
# bcmvn_kidney <- find.pK(sweep.stats_kidney)
# 
# nExp <- round(ncol(rcc.filt) * 0.075)  # expect 5% doublets
# set.seed(555)
# rcc.filt <- doubletFinder(rcc.filt, pN = 0.25, pK = 0.22, nExp = nExp, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
# DF.name = colnames(rcc.filt@meta.data)[grepl("DF.classifications_0.25_0.22_5787", colnames(rcc.filt@meta.data))]
DF.name = colnames(rcc.DF.join@meta.data)[grepl("doublet_finder", colnames(rcc.DF.join@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(rcc.DF.join, group.by = "orig.ident") + NoAxes(),
                   DimPlot(rcc.DF.join, group.by = DF.name) + NoAxes())

VlnPlot(rcc.DF.all, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

rcc.filt = rcc.DF.join[,rcc.DF.join@meta.data$doublet_finder == "Singlet"]
dim(rcc.filt)

p1o <- DimPlot(object = rcc.filt, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2o <- VlnPlot(object = rcc.filt, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1o | p2o

set.seed(555)
rcc.filt <- rcc.filt %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

p1n <- DimPlot(object = rcc.filt, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2n <- VlnPlot(object = rcc.filt, features = "harmony_1",group.by = "orig.ident", pt.size = .1)

# cairo_pdf(paste0(setdir,"rcc_all_harmony_plots.pdf"),width = 18,height = 13)
# print(ggarrange(p1o,p1n,p2o,p2n,nrow=2,ncol=2))
# dev.off()

set.seed(555)
rcc.filt <- rcc.filt %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

# Visualization
orig.stim <- DimPlot(rcc.filt, reduction = "umap", label = TRUE, repel = TRUE,
                     split.by = "orig.ident",ncol = 3,
                     group.by = 'seurat_clusters',pt.size = 0.5,label.size = 5,order = T) + 
  ggtitle("By seurat clusters")

# cairo_pdf(paste0(setdir,"rcc_all_seurat_clusters.pdf"),width = 22,height = 20)
# print(orig.stim)
# dev.off()

condition <- DimPlot(rcc.filt, reduction = "umap",pt.size = 0.5,
                     group.by = "orig.ident",label = T,repel = T,
                     label.size = 5,order=T) + 
  ggtitle("By Tumor Normal")
condition

clusters <- DimPlot(rcc.filt, reduction = "umap",pt.size = 0.5,
                    group.by = "seurat_clusters",label = T,repel = T,
                    label.size = 5,order=T) + 
  ggtitle("By UMAP Clusters")
clusters

# cairo_pdf(paste0(setdir,"rcc_all_clusters.pdf"),width = 14,height = 11)
# print(ggarrange(condition,clusters,nrow=1))
# dev.off()

# Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc.filt) <- "RNA"
rcc.filt <- ScaleData(rcc.filt,verbose = T,features = rownames(rcc.filt))
es.max = sctype_score(scRNAseqData = rcc.filt@assays$RNA$scale.data, 
                      scaled = TRUE, 
                      gs = mrkr.list1)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc.filt@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc.filt@meta.data[rcc.filt@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc.filt@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
rcc.filt$ca9 <- rcc.filt@assays$RNA$scale.data['CA9',]

### UMAP
rcc.filt@meta.data$cellassign1 = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  rcc.filt@meta.data$cellassign1[rcc.filt@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
rcc.filt$cellassign1 <- ifelse(rcc.filt$cellassign1 == 'Proximal tubular cell' & rcc.filt$ca9 > 0,'Proximal Tubular cell + CA9',rcc.filt$cellassign1)

custom1 <- DimPlot(rcc.filt, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",ncol=3,
                   group.by = 'cellassign',pt.size = 0.3,label.size = 5,order = T) + 
  ggtitle("By cell type")
custom1

custom2 <- DimPlot(rcc.filt, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = 'orig.ident', ncol = 3,
                   group.by = 'seurat_clusters',pt.size = 0.3,label.size = 5,order = T) + 
  ggtitle("By cluster")
custom2

cowplot::plot_grid(ncol = 2, DimPlot(rcc.filt, group.by = "orig.ident") + NoAxes(),
                   DimPlot(rcc.filt, group.by = "cellassign") + NoAxes())

# Find Markers
all.markers <- FindAllMarkers(rcc.filt,
                              min.pct = 0.25,
                              min.diff.pct = 0.25)
pos.markers <- FindAllMarkers(rcc.filt,
                              min.pct = 0.25,
                              min.diff.pct = 0.25,
                              verbose = F,
                              only.pos = T)
pos.markers <- order(pos.markers,'cluster','pct.1')
pos.marker.top15 <- pos.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

feat.plots <- FeaturePlot(rcc.filt,
            features = c('RHCG'),
            reduction = 'umap',
            #split.by = 'orig.ident',
            combine = T)
wrap_plots(feat.plots, ncol = 3, nrow = 4)

dim.plots <- DimPlot(rcc.filt,
                     reduction = 'umap')

# cairo_pdf(paste0(setdir,"RCC all plots.pdf"),width = 22,height = 20)
# print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
# dev.off()
# 
# cairo_pdf(paste0(setdir,"RCC all UMAP cell clusters.pdf"),width = 22,height = 15)
# print(custom1)
# dev.off()
# 
# cairo_pdf(paste0(setdir,"RCC all VEGFA CA9.pdf"),width = 30,height = 14)
# print(FeaturePlot(rcc.filt,
#                   repel = T,order = T,
#                   features = c("VEGFA","CA9","EGFR"),
#                   split.by = "orig.ident",
#                   pt.size = 0.3))
# dev.off()

rcc.filt@meta.data$stim.new <- paste0(rcc.filt@meta.data$orig.ident)
rcc.filt@meta.data$stim.new <- gsub("N","",rcc.filt@meta.data$stim.new)
rcc.filt@meta.data$stim.new <- gsub("T1","",rcc.filt@meta.data$stim.new)
rcc.filt@meta.data$stim.new <- gsub("T2","",rcc.filt@meta.data$stim.new)

rcc.filt@meta.data$tissue.type <- factor(rcc.filt@meta.data$orig.ident,
                                         levels = c("RCC1N","RCC1T1","RCC1T2",
                                                    "RCC2N","RCC2T1","RCC2T2",
                                                    "RCC3N","RCC3T1","RCC3T2"),
                                         labels = c(rep(c("normal","tumor","tumor"),3)))

custom3 <- DimPlot(rcc.filt, reduction = "umap", label = TRUE, repel = TRUE,split.by = "tissue.type",
                   group.by = 'cellassign',pt.size = 0.3,label.size = 5,order = T) + 
  ggtitle("By cell type")
custom3

# cairo_pdf(paste0(setdir,"RCC all UMAP cell clusters by tissue type.pdf"),width = 18,height = 9)
# print(custom3)
# dev.off()

FeaturePlot(rcc.filt,
            repel = T,order = T,
            features = c("VEGFA","CA9","EGFR","IGF2BP3"),
            split.by = "tissue.type",
            cols = c("gray","red"),
            pt.size = 0.5)
#####

# Combine all RCCs
#####
rcc.all <- merge(rcc5, y = c(rcc9,rcc10),
                 add.cell.ids = c("RCC5","RCC9","RCC10"))%>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = rcc.all@var.genes, npcs = 20, verbose = FALSE)

allRCC.pre <- DimPlot(rcc.all,
                      reduction = 'pca',
                      group.by = 'orig.ident',
                      label = T,
                      order = T,
                      repel = T) +
  ggtitle('All RCC Pre Harmony PCA')

rcc.all <- RunHarmony(rcc.all,
                      'orig.ident',
                      verbose = F)
rcc.all <- FindNeighbors(rcc.all, dims = 1:20, verbose = F) %>%
  FindClusters(resolution = 0.8, verbose = F) %>%
  RunUMAP(dims = 1:20, verbose = F)

allRCC.post <- DimPlot(rcc.all,
                       reduction = 'harmony',
                       group.by = 'orig.ident',
                       label = T,
                       order = T,
                       repel =T) +
  ggtitle('All RCC Post Harmony PCA')
### Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc.all) <- "RNA" # SCT
rcc.all <- ScaleData(rcc.all,verbose = T,features = rownames(rcc.all))
es.max = sctype_score(scRNAseqData = rcc.all@assays$RNA$scale.data, #SCT
                      scaled = TRUE,
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc.all@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc.all@meta.data[rcc.all@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc.all@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

### UMAP
rcc.all@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  rcc.all@meta.data$cellassign[rcc.all@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

allRCC.tissue <- DimPlot(rcc.all,
                  reduction = 'umap',
                  group.by = 'orig.ident',
                  label = T,
                  order = T,
                  repel = T,
                  cols = c('yellow','orange','red','turquoise','blue','purple','navy','lightblue','green','darkgreen')) +
  ggtitle('All RCC by Tissue Origin')
allRCC.cluster <- DimPlot(rcc.all,
                          reduction = 'umap',
                          group.by = 'seurat_clusters',
                          label = T,
                          repel = T,
                          order = T) +
  ggtitle('All RCC by Cluster')
allRCC.cellassign <- DimPlot(rcc.all,
                             reduction = 'umap',
                             group.by = 'cellassign',
                             label = T,
                             repel = T,
                             order = T) +
  ggtitle('All RCC by Cell Assign')
allRCC.ca9 <- FeaturePlot(rcc.all,
                          features = 'CA9',
                          reduction = 'umap',
                          label = T,
                          repel = T,
                          order = T)
allRCC.cxcl14 <- FeaturePlot(rcc.all,
                             features = 'CXCL14',
                             reduction = 'umap',
                             label = T,
                             repel = T,
                             order = T)
# Derive CCP EMT 15G scores
#####
myr.15 <- c("IGF2BP3","PTK6" ,   "BOLA3" ,  "UQCRH" ,  "TFCP2L1" ,"SLC7A8" , "SLC16A1",
            "PLIN5"  , "TLCD5", "SLC7A1" , "TMEM81",  "RPL22L1", "ELOB"  ,  "IQCC"   ,
            "TAGLN2")
ccp <-  c("FOXM1","ASPM","TK1","PRC1","CDC20",
          "BUB1B",	"PBK",	"DTL",
          "CDKN3","RRM2","ASF1B","CEP55",
          "CDK1","DLGAP5","SKA1","RAD51",
          "KIF11","BIRC5","RAD54L","CENPM",
          "PCLAF","KIF20A","PTTG1","CDCA8",
          "NUSAP1","PLK1","CDCA3","ORC6",
          "CENPF","TOP2A","MCM10")
emt <- c("LRRC15",  "GAS1",    "COL11A1", "PTX3",    "COL7A1",  
         "COL8A2","LOXL2",   "SFRP4",   "PCOLCE",  "MXRA5",   
         "LOXL1",  "DKK1",    "CTHRC1",  "COL6A3",  "IL6",     
         "PRRX1",   "COL1A1",  "FBN2",    "SERPINE1",
         "COMP",    "TFPI2",   "LUM")
expression_matrix <- rcc.filt@assays$RNA@scale.data[c(myr.15,ccp,emt),] %>% as.matrix()
all(colnames(expression_matrix)==rownames(rcc.filt@meta.data))

rcc.filt$CCP.zscore <- as.numeric(rowSums(t(expression_matrix[rownames(expression_matrix) %in% ccp,])))
rcc.filt$EMT.zscore <- as.numeric(rowSums(t(expression_matrix[rownames(expression_matrix) %in% emt,])))
rcc.filt$G15.zscore <- as.numeric(rowSums(t(expression_matrix[rownames(expression_matrix) %in% myr.15,])))

# rcc.filt <- AddModuleScore(rcc.filt,features = list(ccp),name = "CCPscore",nbin = 24,ctrl=100)
# rcc.filt <- AddModuleScore(rcc.filt,features = list(emt),name = "EMTscore",nbin = 24,ctrl=100)
# rcc.filt <- AddModuleScore(rcc.filt,features = list(myr.15),name = "UM15Gscore",nbin = 24,ctrl=100)

FeaturePlot(rcc.filt,
            repel = T,order = T,
            features = c("CCP.zscore","EMT.zscore","G15.zscore"),
            split.by = "tissue.type",
            cols = c("gray","red"),
            pt.size = 0.5)

rcc.filt$Newgroup <- factor(paste0(rcc.filt$tissue.type,"_",rcc.filt$cellassign,sep=""))
dfcomp <- rcc.filt@meta.data

mycomp <- list(c("normal","tumor"))
cairo_pdf("/home/kzayne/ccRCC_ssRNAseq/15Gscores_tumor_normal.pdf",width = 11,height = 9)
dfcomp %>%
  ggplot(aes(x=tissue.type, y=CCP.zscore, fill=tissue.type)) +
  geom_boxplot(alpha=0.8,outlier.shape=NA,width=0.5)+
  ylim(min(dfcomp$CCP.zscore),max(dfcomp$CCP.zscore))+
  scale_fill_manual(values=c("#A3A500","#E76BF3")) +
  scale_y_continuous(expand = expansion(mult = c(0.2,0.35)))+
  # geom_jitter(color="black",size=0.8,alpha=0.9,width =0.25) +
  theme_classic()+theme(text = element_text(size = 16))+
  theme(axis.text.x = element_text(size=10,angle=0))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  ggtitle("") + ylab("15G score")+xlab("")+
  facet_wrap(~cellassign,scales="free")+
  theme(legend.position = "top")+
  stat_compare_means(comparisons=mycomp,method="wilcox",size=5)
dev.off()
#####




# Identify conserved cell type markers
# For performing differential expression after integration, we switch back to the original
# data
#####
DefaultAssay(rcc1.combined.sct) <- "RNA"
rcc1.combined.sct$cellassign <- factor(rcc1.combined.sct$cellassign)
table(rcc1.combined.sct$cellassign,rcc1.combined.sct$seurat_clusters)
cd4t.stim.markers <- FindConservedMarkers(rcc1.combined.sct,
                                          grouping.var = "orig.ident",
                                          ident.1 = 0,
                                          verbose = FALSE)
macro.stim.markers <- FindConservedMarkers(rcc1.combined.sct,
                                           grouping.var = "orig.ident",
                                           ident.1 = 2,
                                           verbose = FALSE)

# # find markers for every cluster compared to all remaining cells, report only the positive ones
# set.seed(123)
# rcc1.markers <- FindAllMarkers(rcc1,only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.75,test.use = "MAST")
# rcc1.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 2, order_by = avg_log2FC) -> sub.markers
# rcc1.markers <- Add_Pct_Diff(marker_dataframe = rcc1.markers)
# rcc1.markers.sub <- rcc1.markers[rcc1.markers$avg_log2FC > 1.5 & rcc1.markers$pct_diff > 0.2,]

# find all markers of cluster 1
cd4t.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "CD4+ T cell", min.pct = 0.25)
natkil.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "Natural killer cell", min.pct = 0.25)
proxi.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "Proximal tubular cell", min.pct = 0.25)
macro.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "Macrophage", min.pct = 0.25)
endo.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "Endothelial cell", min.pct = 0.25)
neut.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "Neutrophil", min.pct = 0.25)
peric.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "Pericyte", min.pct = 0.25)
fibro.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "Fibroblast", min.pct = 0.25)
epith.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "Epithelial cell", min.pct = 0.25)
bcell.markers <- FindMarkers(rcc1.combined.sct,group.by = "cellassign",ident.1 = "B cell", min.pct = 0.25)
imm.markers <- FindMarkers(rcc1,group.by = "customclassif",ident.1 = "Immune cell", min.pct = 0.25)
norm.markers <- FindMarkers(rcc1,group.by = "customclassif",ident.1 = "Normal cell", min.pct = 0.25)

subgenes <- c("LTB","NKG7","GZMA","PTPN22","AIF1","MRPL22","ALDOB","IGFBP7","SPTBN1","MRC1","VEGFA","CD24A","SFRP1","STMN1","VWF")

cairo_pdf("/Users/srnallan/Desktop/Simpa lab/single_cell_ccRCC/RCC1/UMAP cancer markers.pdf",width = 14,height = 9)
my_palette <- colorRampPalette(colors = c("gray","red","orange"))(n = 999)
FeaturePlot(rcc1.aaron,
            repel=F,order = T,
            features = c("VEGFA","CA9","NCAM1","EPCAM"),
            split.by = "orig.ident",
            # group.by="cellassign",
            cols = my_palette,pt.size = 0.5)
dev.off()

# Feature plot CD8A, CD8B, GZMA, GZMB and PRF1
imp.genes <- FeaturePlot(rcc1, features = subgenes,
                         cols = c("gray","red"),pt.size = 0.5,order = F)

imp.genes
# texhaust <- FeaturePlot(rcc1, features = c("Ctla4","Gzma","Gzmb","Lag3","Tigit","Cd274","Tcf1"),
#             cols = c("gray","red"),pt.size = 0.75,order = T)# + ggtitle("TIDE Tcell exhaustion markers")
# texhaust

alex <- FeaturePlot(rcc1, features = c("Hnf1a","Eomes","Pdcd1","Btla","Havcr2"),
                    cols = c("gray","red"),pt.size = 0.75,order = T)# + ggtitle("TIDE Tcell exhaustion markers")
alex

cairo_pdf("/Users/srnallan/Desktop/Simpa lab/single_cell_ccRCC/RCC1/RCC1 all UMAP plots.pdf",width = 18,height = 9)
condition | custom | alex
dev.off()

# Heatmap
rcc1@meta.data$customclassif <- factor(rcc1@meta.data$customclassif,levels = c("CD4+ T cell","CD8+ T cell","Neutrophil",
                                                                               "Proximal tubular cell","Pericyte",
                                                                               "Endothelial cell","Macrophage",
                                                                               "Cancer cell","Epithelial cell","B cell",
                                                                               "Immune cell","Normal cell"))
heat <- DoHeatmap(ScaleData(rcc1,vars.to.regress = "percent.mt",features = subgenes), features = subgenes,
                  group.by = "customclassif") +
  scale_fill_gradientn(colors = c("white","gray","red"))
heat

cairo_pdf("/Users/srnallan/Desktop/Simpa lab/single_cell_ccRCC/RCC1/RCC1 all UMAP plots plus heatmap.pdf",width = 23,height = 12)
g1 <- ggarrange(condition,custom,ncol=2,labels=c("A","B"))
g2 <- ggarrange(g1,heat,ncol=1,nrow=2,labels = c("","C"))
ggarrange(g2,imp.genes,nrow = 1,labels = c("","D"))
dev.off()


# 15G heatmap
myr.15 <- c("IGF2BP3","PTK6" ,   "BOLA3" ,  "UQCRH" ,  "TFCP2L1" ,"SLC7A8" , "SLC16A1",
            "PLIN5"  , "TMEM136", "SLC7A1" , "TMEM81",  "RPL22L1", "ELOB"  ,  "IQCC"   ,
            "TAGLN2")
heat15 <- DoHeatmap(ScaleData(rcc1,vars.to.regress = "percent.mt",features = myr.15), features = myr.15,
                    group.by = "customclassif") +
  scale_fill_gradientn(colors = c("white","gray","red"))
heat15

umap15g <- FeaturePlot(rcc1, features = myr.15,
                       cols = c("gray","red"),pt.size = 0.75,order = T) + # + ggtitle("TIDE Tcell exhaustion markers")+
  patchwork::plot_layout(ncol = 3, nrow = 5)
umap15g
cairo_pdf("/Users/srnallan/Desktop/Simpa lab/single_cell_ccRCC/RCC1/RCC1 15G umap.pdf",width = 20,height = 12)
ggarrange(condition,custom,umap15g,ncol=3,labels=c("A","B","C"),widths = c(1,1,2))
dev.off()

# EMT
emt.sub <- read.csv("/Users/srnallan/Desktop/Simpa lab/IVC Thrombectomy/TCGA/EMT_scoring_UP_060421.csv",header=T)
emt.sub <- as.character(emt.sub$Gene.UP)
umapemt <- FeaturePlot(rcc1, features = emt.sub,
                       cols = c("gray","red"),pt.size = 0.75,order = T)+
  patchwork::plot_layout(ncol = 4, nrow = 6)
umapemt
cairo_pdf("/Users/srnallan/Desktop/Simpa lab/single_cell_ccRCC/RCC1/RCC1 EMT umap.pdf",width = 20,height = 12)
ggarrange(condition,custom,umapemt,ncol=3,labels=c("A","B","C"),widths = c(1,1,2))
dev.off()


# Identify differential expressed genes across conditions
rcc1.combined.sct$celltype.stim <- paste(gsub(" ","_",rcc1.combined.sct$cellassign), rcc1.combined.sct$orig.ident,
                                         sep = "_")
Idents(rcc1.combined.sct) <- "celltype.stim"

rcc1.combined.sct <- PrepSCTFindMarkers(rcc1.combined.sct)
b.interferon.response <- FindMarkers(rcc1.combined.sct, assay = "SCT", ident.1 = "B_STIM", ident.2 = "B_CTRL",
                                     verbose = FALSE)
head(b.interferon.response, n = 15)

library(copykat)
cpykat.1n <- copykat(GetAssayData(object = rcc1.aaron[, WhichCells(rcc1.aaron, ident = c(12,26))], slot = "counts"),
                     id.type = "S", win.size=25, KS.cut=0.1,
                     sam.name="rcc1", distance="euclidean",
                     # norm.cell.names=c(rownames(rcc1@meta.data)[grepl("RCC1N",rownames(rcc1@meta.data))]),
                     n.cores=4, output.seg="FALSE", plot.genes="TRUE", genome="hg20")

pred.test <- data.frame(cpykat.1n$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
CNA.test <- data.frame(cpykat.1n$CNAmat)
# colnames(CNA.test) <- gsub("_.*","",colnames(CNA.test))

CNA.3 <- CNA.test[,4:ncol(CNA.test)]
CNA.3[CNA.3 > 0.4] <- 0.4
CNA.3[CNA.3 < -0.4] <- -0.4
CNA.anno <- cpykat.1n$prediction#[rownames(cpykat.1n$prediction) %in% colnames(CNA.3),]
CNA.anno$cell.names <- gsub("_.*","",CNA.anno$cell.names)
row.anno <- CNA.test[,1:3]
row.anno$chrom <- factor(row.anno$chrom)
chromvar <- c(rep(c("gray","black"),11),"gray")
names(chromvar) <- names(summary(row.anno$chrom))
pheatmap(t(CNA.3),cluster_rows = F,cluster_cols = F,
         annotation_row = CNA.anno,annotation_col = data.frame(chrom=row.anno[,1]),
         annotation_colors = list(chrom=chromvar),
         labels_col = "",labels_row = "",color = colorRampPalette(c("blue","lightblue","white","pink","red"))(100))
#####

#########################################################################
# GSEA and DEG Analyses
#########################################################################
library(EnhancedVolcano)
library(fgsea)
library(GSA)
library(doBy)

# Set cutoffs
fdr <- 0.05
fc <- 2

# Create variable for pseudo-bulk sorting
rcc.filt@meta.data$group <- paste0(rcc.filt@meta.data$orig.ident,"_",rcc.filt@meta.data$cellassign)

# Normalize data; accounts for sequencing depth
rcc.norm <- NormalizeData(
  rcc.filt,
  assay = 'RNA',
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose = TRUE
)

# Scale data; zero centers then scales to bring genes into same range
rcc.norm <- ScaleData(rcc.norm,verbose = T,features = rownames(rcc.norm))

# Join layers
rcc.norm <- JoinLayers(rcc.norm)

# Set metadata column to be used for comparisons
Idents(rcc.norm) <- "orig.ident"

# Hallmark cancer
hm.sym <- GSA.read.gmt("/home/kzayne/Salami-Lab-RCC-Organoid-Project/h.all.v7.4.symbols.gmt.txt")    # download all available genesets genesymbols .gmt file from MSigDB
names(hm.sym$genesets) <- hm.sym$geneset.names
hm.sym <- hm.sym$genesets

# Xcell pathways
xcell.sig <- read.csv("/home/kzayne/Salami-Lab-RCC-Organoid-Project/XCell_pathways.csv",header = T)
xcell <- as.list(as.character(xcell.sig$Genes))
names(xcell) <- xcell.sig$Celltype_Source_ID

for(i in 1:length(xcell)){
  xcell[[i]] <- unlist(strsplit(xcell[[i]], "," ))
}

# Subset data for DEA
rcc.norm.proxTubular <- subset(rcc.norm, subset = cellassign == 'Proximal tubular cell')

# RCC1
#####
# Identify DEGs between T1 and Normal
rcc1.T1vN <- FindMarkers(rcc.norm.proxTubular, ident.1 = "RCC1T1", ident.2 = "RCC1N")

# Create Volcano plot
topdeg1 <- rcc1.T1vN
topdeg1$Comparison <- "RCC1 T1 v N"
topdeg1$Direction <- ifelse(topdeg1$avg_log2FC > 0,'Up','Down')
topdeg1$Significant <- ifelse(abs(topdeg1$avg_log2FC) > 1,'Yes','No')

topdeg1$Direction <- ifelse(topdeg1$Significant == 'No','Not Significant',topdeg1$Direction)

rcc1.volcano.t1vn <- EnhancedVolcano(topdeg1,
                                     lab = rownames(topdeg1),
                                     x = 'avg_log2FC',
                                     y = 'p_val_adj',
                                     # xlim = c(-4,4),
                                     #ylim=c(0,1),
                                     title = paste0(unique(topdeg1$Comparison),
                                                    "\n","nUP=",nrow(topdeg1[topdeg1$Direction == 'Up',]),
                                                    "; nDOWN=",nrow(topdeg1[topdeg1$Direction == 'Down',])),
                                     subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                     pCutoff = fdr,ylab = "FDR adjusted p-value",
                                     xlab = "Higher in Normal <---- Fold Change ----> Higher in T1",
                                     FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                     drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC1 T1 v N volcano.pdf"),width = 22,height = 11)
print(rcc1.volcano.t1vn)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg1#[order(top.table1$logFC,decreasing = T),]
ranks1 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks1) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks1>0,"Up","Down")))

set.seed(123)
fgseaRes1 <- fgseaMultilevel(xcell,ranks1,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot1 <- fgseaRes1[,c("pathway","padj","NES")]

fgseplot1$pathway <- gsub("HALLMARK_","",fgseplot1$pathway)

fgseplot1 <- orderBy(~ -NES + padj,fgseplot1)
fgseplot1$padj <- round(fgseplot1$padj,3)
fgseplot1$Significance <- factor(ifelse(fgseplot1$padj<0.05 & fgseplot1$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot1$padj<0.05 & fgseplot1$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T1","Enriched in Normal","Not significant"))
fgseplot1 <- fgseplot1[fgseplot1$Significance!="Not significant",]
# fgseplot1 <- fgseplot1[!grepl("ESTROGEN",fgseplot1$pathway),]
fgseplot1$pathway <- factor(fgseplot1$pathway,levels = c(fgseplot1$pathway))
fill <- c("red","blue")

fg1 <- ggplot(data = fgseplot1, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in N <-- NES --> Enriched in T1")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg1$Comparison),",nGenes=",nrow(topdeg1))) +
  scale_x_discrete(limits = rev(levels(fgseplot1$pathway)))

cairo_pdf(paste0(setdir,"RCC1 T1 v N fgsea.pdf"),width = 22,height = 11)
print(fg1)
dev.off()

# Identify DEGs between T2 and Normal
rcc1.T2vN <- FindMarkers(rcc.norm, ident.1 = "RCC1T2", ident.2 = "RCC1N")

# Create Volcano plot
topdeg1 <- rcc1.T2vN
topdeg1$Comparison <- "RCC1 T2 v N"
topdeg1$Direction <- ifelse(topdeg1$avg_log2FC > 0,'Up','Down')
topdeg1$Significant <- ifelse(abs(topdeg1$avg_log2FC) > 1,'Yes','No')

topdeg1$Direction <- ifelse(topdeg1$Significant == 'No','Not Significant',topdeg1$Direction)

rcc1.volcano.t2vn <- EnhancedVolcano(topdeg1,
                                     lab = rownames(topdeg1),
                                     x = 'avg_log2FC',
                                     y = 'p_val_adj',
                                     # xlim = c(-4,4),
                                     #ylim=c(0,1),
                                     title = paste0(unique(topdeg1$Comparison),
                                                    "\n","nUP=",nrow(topdeg1[topdeg1$Direction == 'Up',]),
                                                    "; nDOWN=",nrow(topdeg1[topdeg1$Direction == 'Down',])),
                                     subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                     pCutoff = fdr,ylab = "FDR adjusted p-value",
                                     xlab = "Higher in Normal <---- Fold Change ----> Higher in T2",
                                     FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                     drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC1 T2 v N volcano.pdf"),width = 22,height = 11)
print(rcc1.volcano.t2vn)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg1#[order(top.table1$logFC,decreasing = T),]
ranks1 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks1) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks1>0,"Up","Down")))

set.seed(123)
fgseaRes2 <- fgseaMultilevel(xcell,ranks1,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot2 <- fgseaRes2[,c("pathway","padj","NES")]

fgseplot2$pathway <- gsub("HALLMARK_","",fgseplot2$pathway)

fgseplot2 <- orderBy(~ -NES + padj,fgseplot2)
fgseplot2$padj <- round(fgseplot2$padj,3)
fgseplot2$Significance <- factor(ifelse(fgseplot2$padj<0.05 & fgseplot2$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot2$padj<0.05 & fgseplot2$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T2","Enriched in Normal","Not significant"))
fgseplot2 <- fgseplot2[fgseplot2$Significance!="Not significant",]
# fgseplot2 <- fgseplot2[!grepl("ESTROGEN",fgseplot2$pathway),]
fgseplot2$pathway <- factor(fgseplot2$pathway,levels = c(fgseplot2$pathway))
fill <- c("red","blue")

fg1 <- ggplot(data = fgseplot2, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in N <-- NES --> Enriched in T2")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg1$Comparison),",nGenes=",nrow(topdeg1))) +
  scale_x_discrete(limits = rev(levels(fgseplot2$pathway)))

cairo_pdf(paste0(setdir,"RCC1 T2 v N fgsea.pdf"),width = 22,height = 11)
print(fg1)
dev.off()

# Identify DEGs between T2 and Normal
rcc1.T2vT1 <- FindMarkers(rcc.norm, ident.1 = "RCC1T2", ident.2 = "RCC1T1")

# Create Volcano plot
topdeg1 <- rcc1.T2vT1
topdeg1$Comparison <- "RCC1 T2 v T1"
topdeg1$Direction <- ifelse(topdeg1$avg_log2FC > 0,'Up','Down')
topdeg1$Significant <- ifelse(abs(topdeg1$avg_log2FC) > 1,'Yes','No')

topdeg1$Direction <- ifelse(topdeg1$Significant == 'No','Not Significant',topdeg1$Direction)

rcc1.volcano.t2vt1 <- EnhancedVolcano(topdeg1,
                                      lab = rownames(topdeg1),
                                      x = 'avg_log2FC',
                                      y = 'p_val_adj',
                                      # xlim = c(-4,4),
                                      #ylim=c(0,1),
                                      title = paste0(unique(topdeg1$Comparison),
                                                     "\n","nUP=",nrow(topdeg1[topdeg1$Direction == 'Up',]),
                                                     "; nDOWN=",nrow(topdeg1[topdeg1$Direction == 'Down',])),
                                      subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                      pCutoff = fdr,ylab = "FDR adjusted p-value",
                                      xlab = "Higher in T1 <---- Fold Change ----> Higher in T2",
                                      FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                      drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC1 T2 v T1 volcano.pdf"),width = 22,height = 11)
print(rcc1.volcano.t2vt1)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg1#[order(top.table1$logFC,decreasing = T),]
ranks1 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks1) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks1>0,"Up","Down")))

set.seed(123)
fgseaRes3 <- fgseaMultilevel(xcell,ranks1,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot3 <- fgseaRes3[,c("pathway","padj","NES")]

fgseplot3$pathway <- gsub("HALLMARK_","",fgseplot3$pathway)

fgseplot3 <- orderBy(~ -NES + padj,fgseplot3)
fgseplot3$padj <- round(fgseplot3$padj,3)
fgseplot3$Significance <- factor(ifelse(fgseplot3$padj<0.05 & fgseplot3$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot3$padj<0.05 & fgseplot3$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T2","Enriched in T1","Not significant"))
fgseplot3 <- fgseplot3[fgseplot3$Significance!="Not significant",]
# fgseplot3 <- fgseplot3[!grepl("ESTROGEN",fgseplot3$pathway),]
fgseplot3$pathway <- factor(fgseplot3$pathway,levels = c(fgseplot3$pathway))
fill <- c("red","blue")

fg1 <- ggplot(data = fgseplot3, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in T1 <-- NES --> Enriched in T2")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg1$Comparison),",nGenes=",nrow(topdeg1))) +
  scale_x_discrete(limits = rev(levels(fgseplot3$pathway)))

cairo_pdf(paste0(setdir,"RCC1 T2 v T1 fgsea.pdf"),width = 22,height = 11)
print(fg1)
dev.off()
#####

# RCC2
#####
# Identify DEGs between T1 and Normal
rcc2.T1vN <- FindMarkers(rcc.norm, ident.1 = "RCC2T1", ident.2 = "RCC2N")

# Create Volcano plot
topdeg2 <- rcc2.T1vN
topdeg2$Comparison <- "RCC2 T1 v N"
topdeg2$Direction <- ifelse(topdeg2$avg_log2FC > 0,'Up','Down')
topdeg2$Significant <- ifelse(abs(topdeg2$avg_log2FC) > 1,'Yes','No')

topdeg2$Direction <- ifelse(topdeg2$Significant == 'No','Not Significant',topdeg2$Direction)

rcc2.volcano.t1vn <- EnhancedVolcano(topdeg2,
                                     lab = rownames(topdeg2),
                                     x = 'avg_log2FC',
                                     y = 'p_val_adj',
                                     # xlim = c(-4,4),
                                     #ylim=c(0,1),
                                     title = paste0(unique(topdeg2$Comparison),
                                                    "\n","nUP=",nrow(topdeg2[topdeg2$Direction == 'Up',]),
                                                    "; nDOWN=",nrow(topdeg2[topdeg2$Direction == 'Down',])),
                                     subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                     pCutoff = fdr,ylab = "FDR adjusted p-value",
                                     xlab = "Higher in Normal <---- Fold Change ----> Higher in T1",
                                     FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                     drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC2 T1 v N volcano.pdf"),width = 22,height = 11)
print(rcc2.volcano.t1vn)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg2#[order(top.table1$logFC,decreasing = T),]
ranks2 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks2) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks2>0,"Up","Down")))

set.seed(123)
fgseaRes4 <- fgseaMultilevel(xcell,ranks2,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot4 <- fgseaRes4[,c("pathway","padj","NES")]

fgseplot4$pathway <- gsub("HALLMARK_","",fgseplot4$pathway)

fgseplot4 <- orderBy(~ -NES + padj,fgseplot4)
fgseplot4$padj <- round(fgseplot4$padj,3)
fgseplot4$Significance <- factor(ifelse(fgseplot4$padj<0.05 & fgseplot4$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot4$padj<0.05 & fgseplot4$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T1","Enriched in Normal","Not significant"))
fgseplot4 <- fgseplot4[fgseplot4$Significance!="Not significant",]
# fgseplot4 <- fgseplot4[!grepl("ESTROGEN",fgseplot4$pathway),]
fgseplot4$pathway <- factor(fgseplot4$pathway,levels = c(fgseplot4$pathway))
fill <- c("red","blue")

fg2 <- ggplot(data = fgseplot4, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in N <-- NES --> Enriched in T1")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg2$Comparison),",nGenes=",nrow(topdeg2))) +
  scale_x_discrete(limits = rev(levels(fgseplot4$pathway)))

cairo_pdf(paste0(setdir,"RCC2 T1 v N fgsea.pdf"),width = 22,height = 11)
print(fg2)
dev.off()

# Identify DEGs between T2 and Normal
rcc2.T2vN <- FindMarkers(rcc.norm, ident.1 = "RCC2T2", ident.2 = "RCC2N")

# Create Volcano plot
topdeg2 <- rcc2.T2vN
topdeg2$Comparison <- "RCC2 T2 v N"
topdeg2$Direction <- ifelse(topdeg2$avg_log2FC > 0,'Up','Down')
topdeg2$Significant <- ifelse(abs(topdeg2$avg_log2FC) > 1,'Yes','No')

topdeg2$Direction <- ifelse(topdeg2$Significant == 'No','Not Significant',topdeg2$Direction)

rcc2.volcano.t2vn <- EnhancedVolcano(topdeg2,
                                     lab = rownames(topdeg2),
                                     x = 'avg_log2FC',
                                     y = 'p_val_adj',
                                     # xlim = c(-4,4),
                                     #ylim=c(0,1),
                                     title = paste0(unique(topdeg2$Comparison),
                                                    "\n","nUP=",nrow(topdeg2[topdeg2$Direction == 'Up',]),
                                                    "; nDOWN=",nrow(topdeg2[topdeg2$Direction == 'Down',])),
                                     subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                     pCutoff = fdr,ylab = "FDR adjusted p-value",
                                     xlab = "Higher in Normal <---- Fold Change ----> Higher in T2",
                                     FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                     drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC2 T2 v N volcano.pdf"),width = 22,height = 11)
print(rcc2.volcano.t2vn)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg2#[order(top.table1$logFC,decreasing = T),]
ranks2 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks2) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks2>0,"Up","Down")))

set.seed(123)
fgseaRes5 <- fgseaMultilevel(xcell,ranks2,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot5 <- fgseaRes5[,c("pathway","padj","NES")]

fgseplot5$pathway <- gsub("HALLMARK_","",fgseplot5$pathway)

fgseplot5 <- orderBy(~ -NES + padj,fgseplot5)
fgseplot5$padj <- round(fgseplot5$padj,3)
fgseplot5$Significance <- factor(ifelse(fgseplot5$padj<0.05 & fgseplot5$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot5$padj<0.05 & fgseplot5$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T2","Enriched in Normal","Not significant"))
fgseplot5 <- fgseplot5[fgseplot5$Significance!="Not significant",]
# fgseplot5 <- fgseplot5[!grepl("ESTROGEN",fgseplot5$pathway),]
fgseplot5$pathway <- factor(fgseplot5$pathway,levels = c(fgseplot5$pathway))
fill <- c("red","blue")

fg2 <- ggplot(data = fgseplot5, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in N <-- NES --> Enriched in T2")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg2$Comparison),",nGenes=",nrow(topdeg2))) +
  scale_x_discrete(limits = rev(levels(fgseplot5$pathway)))

cairo_pdf(paste0(setdir,"RCC2 T2 v N fgsea.pdf"),width = 22,height = 11)
print(fg2)
dev.off()

# Identify DEGs between T2 and Normal
rcc2.T2vT1 <- FindMarkers(rcc.norm, ident.1 = "RCC2T2", ident.2 = "RCC2T1")

# Create Volcano plot
topdeg2 <- rcc2.T2vT1
topdeg2$Comparison <- "RCC2 T2 v T1"
topdeg2$Direction <- ifelse(topdeg2$avg_log2FC > 0,'Up','Down')
topdeg2$Significant <- ifelse(abs(topdeg2$avg_log2FC) > 1,'Yes','No')

topdeg2$Direction <- ifelse(topdeg2$Significant == 'No','Not Significant',topdeg2$Direction)

rcc2.volcano.t2vt1 <- EnhancedVolcano(topdeg2,
                                      lab = rownames(topdeg2),
                                      x = 'avg_log2FC',
                                      y = 'p_val_adj',
                                      # xlim = c(-4,4),
                                      #ylim=c(0,1),
                                      title = paste0(unique(topdeg2$Comparison),
                                                     "\n","nUP=",nrow(topdeg2[topdeg2$Direction == 'Up',]),
                                                     "; nDOWN=",nrow(topdeg2[topdeg2$Direction == 'Down',])),
                                      subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                      pCutoff = fdr,ylab = "FDR adjusted p-value",
                                      xlab = "Higher in T1 <---- Fold Change ----> Higher in T2",
                                      FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                      drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC2 T2 v T1 volcano.pdf"),width = 22,height = 11)
print(rcc2.volcano.t2vt1)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg2#[order(top.table1$logFC,decreasing = T),]
ranks2 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks2) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks2>0,"Up","Down")))

set.seed(123)
fgseaRes6 <- fgseaMultilevel(xcell,ranks2,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot6 <- fgseaRes6[,c("pathway","padj","NES")]

fgseplot6$pathway <- gsub("HALLMARK_","",fgseplot6$pathway)

fgseplot6 <- orderBy(~ -NES + padj,fgseplot6)
fgseplot6$padj <- round(fgseplot6$padj,3)
fgseplot6$Significance <- factor(ifelse(fgseplot6$padj<0.05 & fgseplot6$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot6$padj<0.05 & fgseplot6$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T2","Enriched in T1","Not significant"))
fgseplot6 <- fgseplot6[fgseplot6$Significance!="Not significant",]
# fgseplot6 <- fgseplot6[!grepl("ESTROGEN",fgseplot6$pathway),]
fgseplot6$pathway <- factor(fgseplot6$pathway,levels = c(fgseplot6$pathway))
fill <- c("red","blue")

fg2 <- ggplot(data = fgseplot6, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in T1 <-- NES --> Enriched in T2")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg2$Comparison),",nGenes=",nrow(topdeg2))) +
  scale_x_discrete(limits = rev(levels(fgseplot6$pathway)))

cairo_pdf(paste0(setdir,"RCC2 T2 v T1 fgsea.pdf"),width = 22,height = 11)
print(fg2)
dev.off()
#####

# RCC3
#####
# Identify DEGs between T1 and Normal
rcc3.T1vN <- FindMarkers(rcc.norm, ident.1 = "RCC3T1", ident.2 = "RCC3N")

# Create Volcano plot
topdeg3 <- rcc3.T1vN
topdeg3$Comparison <- "RCC3 T1 v N"
topdeg3$Direction <- ifelse(topdeg3$avg_log2FC > 0,'Up','Down')
topdeg3$Significant <- ifelse(abs(topdeg3$avg_log2FC) > 1,'Yes','No')

topdeg3$Direction <- ifelse(topdeg3$Significant == 'No','Not Significant',topdeg3$Direction)

rcc3.volcano.t1vn <- EnhancedVolcano(topdeg3,
                                     lab = rownames(topdeg3),
                                     x = 'avg_log2FC',
                                     y = 'p_val_adj',
                                     # xlim = c(-4,4),
                                     #ylim=c(0,1),
                                     title = paste0(unique(topdeg3$Comparison),
                                                    "\n","nUP=",nrow(topdeg3[topdeg3$Direction == 'Up',]),
                                                    "; nDOWN=",nrow(topdeg3[topdeg3$Direction == 'Down',])),
                                     subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                     pCutoff = fdr,ylab = "FDR adjusted p-value",
                                     xlab = "Higher in Normal <---- Fold Change ----> Higher in T1",
                                     FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                     drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC3 T1 v N volcano.pdf"),width = 22,height = 11)
print(rcc3.volcano.t1vn)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg3#[order(top.table1$logFC,decreasing = T),]
ranks3 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks3) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks3>0,"Up","Down")))

set.seed(123)
fgseaRes7 <- fgseaMultilevel(xcell,ranks3,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot7 <- fgseaRes7[,c("pathway","padj","NES")]

fgseplot7$pathway <- gsub("HALLMARK_","",fgseplot7$pathway)

fgseplot7 <- orderBy(~ -NES + padj,fgseplot7)
fgseplot7$padj <- round(fgseplot7$padj,3)
fgseplot7$Significance <- factor(ifelse(fgseplot7$padj<0.05 & fgseplot7$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot7$padj<0.05 & fgseplot7$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T1","Enriched in Normal","Not significant"))
fgseplot7 <- fgseplot7[fgseplot7$Significance!="Not significant",]
# fgseplot7 <- fgseplot7[!grepl("ESTROGEN",fgseplot7$pathway),]
fgseplot7$pathway <- factor(fgseplot7$pathway,levels = c(fgseplot7$pathway))
fill <- c("red","blue")

fg3 <- ggplot(data = fgseplot7, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in N <-- NES --> Enriched in T1")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg3$Comparison),",nGenes=",nrow(topdeg3))) +
  scale_x_discrete(limits = rev(levels(fgseplot7$pathway)))

cairo_pdf(paste0(setdir,"RCC3 T1 v N fgsea.pdf"),width = 22,height = 11)
print(fg3)
dev.off()

# Identify DEGs between T2 and Normal
rcc3.T2vN <- FindMarkers(rcc.norm, ident.1 = "RCC3T2", ident.2 = "RCC3N")

# Create Volcano plot
topdeg3 <- rcc3.T2vN
topdeg3$Comparison <- "RCC3 T2 v N"
topdeg3$Direction <- ifelse(topdeg3$avg_log2FC > 0,'Up','Down')
topdeg3$Significant <- ifelse(abs(topdeg3$avg_log2FC) > 1,'Yes','No')

topdeg3$Direction <- ifelse(topdeg3$Significant == 'No','Not Significant',topdeg3$Direction)

rcc3.volcano.t2vn <- EnhancedVolcano(topdeg3,
                                     lab = rownames(topdeg3),
                                     x = 'avg_log2FC',
                                     y = 'p_val_adj',
                                     # xlim = c(-4,4),
                                     #ylim=c(0,1),
                                     title = paste0(unique(topdeg3$Comparison),
                                                    "\n","nUP=",nrow(topdeg3[topdeg3$Direction == 'Up',]),
                                                    "; nDOWN=",nrow(topdeg3[topdeg3$Direction == 'Down',])),
                                     subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                     pCutoff = fdr,ylab = "FDR adjusted p-value",
                                     xlab = "Higher in Normal <---- Fold Change ----> Higher in T2",
                                     FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                     drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC3 T2 v N volcano.pdf"),width = 22,height = 11)
print(rcc3.volcano.t2vn)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg3#[order(top.table1$logFC,decreasing = T),]
ranks3 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks3) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks3>0,"Up","Down")))

set.seed(123)
fgseaRes8 <- fgseaMultilevel(xcell,ranks3,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot8 <- fgseaRes8[,c("pathway","padj","NES")]

fgseplot8$pathway <- gsub("HALLMARK_","",fgseplot8$pathway)

fgseplot8 <- orderBy(~ -NES + padj,fgseplot8)
fgseplot8$padj <- round(fgseplot8$padj,3)
fgseplot8$Significance <- factor(ifelse(fgseplot8$padj<0.05 & fgseplot8$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot8$padj<0.05 & fgseplot8$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T2","Enriched in Normal","Not significant"))
fgseplot8 <- fgseplot8[fgseplot8$Significance!="Not significant",]
# fgseplot8 <- fgseplot8[!grepl("ESTROGEN",fgseplot8$pathway),]
fgseplot8$pathway <- factor(fgseplot8$pathway,levels = c(fgseplot8$pathway))
fill <- c("red","blue")

fg3 <- ggplot(data = fgseplot8, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in N <-- NES --> Enriched in T2")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg3$Comparison),",nGenes=",nrow(topdeg3))) +
  scale_x_discrete(limits = rev(levels(fgseplot8$pathway)))

cairo_pdf(paste0(setdir,"RCC3 T2 v N fgsea.pdf"),width = 22,height = 11)
print(fg3)
dev.off()

# Identify DEGs between T2 and Normal
rcc3.T2vT1 <- FindMarkers(rcc.norm, ident.1 = "RCC3T2", ident.2 = "RCC3T1")

# Create Volcano plot
topdeg3 <- rcc3.T2vT1
topdeg3$Comparison <- "RCC3 T2 v T1"
topdeg3$Direction <- ifelse(topdeg3$avg_log2FC > 0,'Up','Down')
topdeg3$Significant <- ifelse(abs(topdeg3$avg_log2FC) > 1,'Yes','No')

topdeg3$Direction <- ifelse(topdeg3$Significant == 'No','Not Significant',topdeg3$Direction)

rcc3.volcano.t2vt1 <- EnhancedVolcano(topdeg3,
                                      lab = rownames(topdeg3),
                                      x = 'avg_log2FC',
                                      y = 'p_val_adj',
                                      # xlim = c(-4,4),
                                      #ylim=c(0,1),
                                      title = paste0(unique(topdeg3$Comparison),
                                                     "\n","nUP=",nrow(topdeg3[topdeg3$Direction == 'Up',]),
                                                     "; nDOWN=",nrow(topdeg3[topdeg3$Direction == 'Down',])),
                                      subtitle = paste0("Limma; FDR<",fdr,"; absolute LogFC >",fc),
                                      pCutoff = fdr,ylab = "FDR adjusted p-value",
                                      xlab = "Higher in T1 <---- Fold Change ----> Higher in T2",
                                      FCcutoff = fc,pointSize = 4,labSize = 4,maxoverlapsConnectors = 30,
                                      drawConnectors = T,widthConnectors = 0.1,colConnectors = 'grey30')

cairo_pdf(paste0(setdir,"RCC3 T2 v T1 volcano.pdf"),width = 22,height = 11)
print(rcc3.volcano.t2vt1)
dev.off()

# Create Pathway Enrichment plot
# Prep data for pathway enrichment analysis
dfgsea <- topdeg3#[order(top.table1$logFC,decreasing = T),]
ranks3 <- dfgsea[order(-dfgsea[,"avg_log2FC"]),]$avg_log2FC
names(ranks3) <- rownames(dfgsea[order(-dfgsea[,"avg_log2FC"]),])
table(factor(ifelse(ranks3>0,"Up","Down")))

set.seed(123)
fgseaRes9 <- fgseaMultilevel(xcell,ranks3,minSize=10,maxSize = 500,gseaParam = 1,eps=0,nPermSimple = 1000)
fgseplot9 <- fgseaRes9[,c("pathway","padj","NES")]

fgseplot9$pathway <- gsub("HALLMARK_","",fgseplot9$pathway)

fgseplot9 <- orderBy(~ -NES + padj,fgseplot9)
fgseplot9$padj <- round(fgseplot9$padj,3)
fgseplot9$Significance <- factor(ifelse(fgseplot9$padj<0.05 & fgseplot9$NES>0,"Up; FDR < 5%",
                                        ifelse(fgseplot9$padj<0.05 & fgseplot9$NES<0,"Down; FDR < 5%","Not significant")),
                                 levels = c("Up; FDR < 5%","Down; FDR < 5%","Not significant"),
                                 labels=c("Enriched in T2","Enriched in T1","Not significant"))
fgseplot9 <- fgseplot9[fgseplot9$Significance!="Not significant",]
# fgseplot9 <- fgseplot9[!grepl("ESTROGEN",fgseplot9$pathway),]
fgseplot9$pathway <- factor(fgseplot9$pathway,levels = c(fgseplot9$pathway))
fill <- c("red","blue")

fg3 <- ggplot(data = fgseplot9, aes(x = pathway, y = NES,fill=Significance)) +
  geom_bar(stat="identity",width=.5) + coord_flip()+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  scale_fill_manual(values =fill) + xlab("")+ylab("Enriched in T1 <-- NES --> Enriched in T2")+
  theme(axis.title.x = element_text(size=14))+ theme(legend.text = element_text(size = 15))+
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text.x = element_text(size =  15,colour = "black")) +
  theme(plot.title = element_text(size =  15,colour = "black",face = "bold")) +
  ggtitle(paste0("Hallmark; ", unique(topdeg3$Comparison),",nGenes=",nrow(topdeg3))) +
  scale_x_discrete(limits = rev(levels(fgseplot9$pathway)))

cairo_pdf(paste0(setdir,"RCC3 T2 v T1 fgsea.pdf"),width = 22,height = 11)
print(fg3)
dev.off()
#####

############## Fgseplots combined - bubbleplot
# add comparison column to identify where each fgseaResX came from
fgseaRes1$Comparison <- 'RCC1 T1vN'
fgseaRes2$Comparison <- 'RCC1 T2vN'
fgseaRes3$Comparison <- 'RCC1 T2vT1'
fgseaRes4$Comparison <- 'RCC2 T1vN'
fgseaRes5$Comparison <- 'RCC2 T2vN'
fgseaRes6$Comparison <- 'RCC2 T2vT1'
fgseaRes7$Comparison <- 'RCC3 T1vN'
fgseaRes8$Comparison <- 'RCC3 T2vN'
fgseaRes9$Comparison <- 'RCC3 T2vT1'

tops <- rbind(fgseaRes1,fgseaRes2,fgseaRes3,fgseaRes4,fgseaRes5,
              fgseaRes6,fgseaRes7,fgseaRes8,fgseaRes9)
tops <- tops[tops$padj < 0.05,]
tops <- tops[abs(tops$NES) > 1,]
tops <- tops[,c(1:7,9)]
tops <- tops[grepl("HALLMARK",tops$pathway),]
tops$pathway <- gsub("HALLMARK_","",tops$pathway)
tops$pathway <- gsub("_"," ",tops$pathway)
# tops <- tops[!grepl("ESTROGEN",tops$pathway),]
tops <- tops[order(tops$padj,decreasing = F),]
tops$pathway <- factor(tops$pathway,levels = unique(c(tops$pathway)))

tops$Comparison <- factor(tops$Comparison, levels = c("RCC1 T1vN",
                                                      "RCC1 T2vN",
                                                      "RCC1 T2vT1",
                                                      "RCC2 T1vN",
                                                      "RCC2 T2vN",
                                                      "RCC2 T2vT1",
                                                      "RCC3 T1vN",
                                                      "RCC3 T2vN",
                                                      "RCC3 T2vT1"),
                          labels = c("RCC1 T1vN",
                                     "RCC1 T2vN",
                                     "RCC1 T2vT1",
                                     "RCC2 T1vN",
                                     "RCC2 T2vN",
                                     "RCC2 T2vT1",
                                     "RCC3 T1vN",
                                     "RCC3 T2vN",
                                     "RCC3 T2vT1"))

tops$NES <- as.numeric(round(tops$NES,digits=2))

tops$Direction <- factor(ifelse(tops$padj < 0.05 & tops$NES>0,"High in T1/T2",
                                ifelse(tops$padj < 0.05 & tops$NES<0,"High in N/T1","Not significant")),
                         levels = c("High in T1/T2", "High in N/T1","Not significant"))
tops <- tops[!tops$Direction=="Not significant",]
tops$absNES <- abs(tops$NES)

# ggplot combined fgsea bubble plot
fill <- c("red","#3182bd")

g2 <- ggplot(tops, aes(x = Comparison, y = pathway, size = absNES, fill = Direction)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_manual(values =fill) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text=element_text(size=10,colour = "black")) + scale_y_discrete(limits = rev(levels(factor(tops$pathway)))) +
  scale_size(range = c(3,8),breaks = c(1,1.5,2)) +
  ggtitle("Single Cell RNASeq Xcell") +
  xlab("")+ylab("")+theme_bw() +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 30, vjust = 0, hjust = 0,face="bold",size=20,colour = "black")) +
  theme(axis.text.y = element_text(size =  15,colour = "black")) +
  theme(axis.text=element_text(size=20), legend.text=element_text(size=15))+
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

cairo_pdf(paste0(setdir,"Combined RCC fgsea bubble plot xcell.pdf"),width = 22,height = 11)
print(g2)
dev.off()

################################################################################
# Single Cell Pathway Enrichment - UCell
################################################################################
library(scGSVA)
library(GSA)
library(UCell)
library(Cairo)

# Subset data into individual patient
rcc1.norm <- subset(rcc.norm, subset = stim.new == 'RCC1')
rcc2.norm <- subset(rcc.norm, subset = stim.new == 'RCC2')
rcc3.norm <- subset(rcc.norm, subset = stim.new == 'RCC3')

# Join Layers before running UCell
rcc1.norm.join <- JoinLayers(rcc1.norm)

# Load hm.sym hallmark pathways
hm.sym <- GSA.read.gmt(paste0("h.all.v7.4.symbols.gmt.txt"))    # download all available genesets genesymbols .gmt file from MSigDB
names(hm.sym$genesets) <- hm.sym$geneset.names
hm.sym <- hm.sym$genesets
names(hm.sym) <- gsub("HALLMARK_","",names(hm.sym))
hm.sym <- hm.sym[order(names(hm.sym))]

# Run UCELL
rcc1.norm.U <- AddModuleScore_UCell(rcc1.norm.join,
                                    features = hm.sym)

# Rename names for hm.sym
colnames(rcc1.norm.U@meta.data) <- gsub("_UCell","",colnames(rcc1.norm.U@meta.data))
colnames(rcc1.norm.U@meta.data) <- gsub("HALLMARK_","",colnames(rcc1.norm.U@meta.data))

# Factor
rcc1.norm.U@meta.data$orig.ident <- factor(rcc1.norm.U@meta.data$orig.ident)

# Plot UMAP results
signatures <- c("E2F_TARGETS","G2M_CHECKPOINT","EPITHELIAL_MESENCHYMAL_TRANSITION")
signatures1 <- c("INTERFERON_GAMMA_RESPONSE","INFLAMMATORY_RESPONSE","TNFA_SIGNALING_VIA_NFKB")
fp <- FeaturePlot(rcc1.norm.U,
            reduction = "umap",
            features = "INTERFERON_GAMMA_RESPONSE", #c(signatures1,'ca9'),
            split.by = "orig.ident",
            #cols = c('gray','gray','red'),
            repel = T,
            order = T)
fp <- fp & scale_colour_gradientn(colors = colorRampPalette(c("blue","lightgreen","khaki1","orange","red4"))(5))#c("#2b83ba","#abdda4","#ffffbf","#fc8d59","#d7191c"))(5)#c('gray','lightpink','maroon')) #c('gray', 'red'))

dp <- DimPlot(rcc1.norm.U,
        reduction = 'umap',
        group.by = 'cellassign',
        split.by = 'orig.ident',
        order = T,
        label = T,
        repel = T)
dp <- dp  + theme(legend.position = 'none')#+ xlim(-15,10) + ylim(-10,10)

# Make loop to print pathways
CairoPDF("RCC1 All Pathways.pdf", width = 19, height = 11, onefile = T)
for( i in names(hm.sym)){
  fp <- FeaturePlot(rcc1.norm.U,
                    reduction = "umap",
                    features = i, #c(signatures1,'ca9'),
                    split.by = "orig.ident",
                    #cols = c('gray','gray','red'),
                    repel = T,
                    order = T)
  fp <- fp & scale_colour_gradientn(colors = colorRampPalette(c("blue","lightgreen","khaki1","orange","red4"))(5))#c("#2b83ba","#abdda4","#ffffbf","#fc8d59","#d7191c"))(5))#c('gray','lightpink','maroon')) #c('gray', 'red'))
  
  print(fp / dp)
}
dev.off()
#####
#gsva.anno <- buildAnnot(species = 'human',keytype='SYMBOL',anntype='GO')
# res <- scgsva(rcc1.norm.join, 
#               annot = hm.sym,
#               method = 'UMAP')


# all(rownames(counts(filtered.full))==rownames(y.full$genes))
# e.gsva <- cpm(counts(filtered.full),log = T)
# e.gsva <- e.gsva[,colnames(e.gsva) %in% rownames(df.des2)]
# rmgns <- y.full$genes
# rmgns <- rmgns[!duplicated(rmgns$gene_name),]
# e.gsva <- e.gsva[rownames(e.gsva) %in% rownames(rmgns),]
# all(rownames(e.gsva)==rownames(rmgns))
# rownames(e.gsva) <- rmgns$gene_name
# 
# 
# set.seed(123)
# fg2 <- gsva(e.gsva,all.sig1)
# 
# dfplot <- annotations
# all(rownames(dfplot)==rownames(annotations2))
# fg2 <- fg2[,match(rownames(annotations2),colnames(fg2))]
# dfplot <- dfplot[match(rownames(annotations2),rownames(dfplot)),]
# all(rownames(dfplot)==colnames(fg2))
# rownames(fg2) <- gsub("HALLMARK_","",rownames(fg2))
# rownames(fg2) <- gsub("_"," ",rownames(fg2),fixed = T)
# fg2 <- fg2[rownames(fg2) %in% fgseplot1$pathway,]
# fg2 <- fg2[match(fgseplot1$pathway,rownames(fg2)),colnames(fg2) %in% rownames(df.des2)]
# 
# CairoPDF("/Users/srnallan/Downloads/PRI_MET_revision/ssgsva.pdf",width = 15,height = 8,onefile = T)
# ComplexHeatmap::pheatmap(fg2,color = colorRampPalette(c("#3182bd","white","red"))(100),
#                          annotation_colors = ann_colors,
#                          fontsize_col = 10,cluster_cols = T,cluster_rows = F,
#                          cellheight = 10,cellwidth = 10,
#                          fontsize_row = 12,fontsize = 9,
#                          annotation_col = annotations[rownames(annotations) %in% rownames(df.des2),],
#                          clustering_distance_cols = "correlation",
#                          labels_col = gsub("Myriad","",annotations2[rownames(annotations2) %in% rownames(df.des2),]$Myriad))
# 
# dev.off()
#####

# Heat map of pathways by cell type
#####
library(pheatmap)

# Get matrix of pathway values from Seurat metadata
df.heatmap <- rcc1.norm.U@meta.data[,c(1,7,9:58)]

# Combine scores from each RCC by cellassign type
df.heatmap.ag <- aggregate(df.heatmap, 
                           by = .~orig.ident+cellassign,
                           FUN = mean)
# # Add rownames
 rownames(df.heatmap.ag) <- paste0(df.heatmap.ag$orig.ident,"_",df.heatmap.ag$cellassign)
# rownames(df.heatmap.ag) = NULL

# Sort out annotation columns i.e. orig.ident and cell assign
df.heatmap.anno <- df.heatmap.ag[,c(1,2)]
df.heatmap.ag <- df.heatmap.ag[,-c(1,2)]
df.heatmap.ag <- as.data.frame(t(df.heatmap.ag))

# Filter out low expressed pathways
df.heatmap.ag$sums <- rowSums(df.heatmap.ag)
df.heatmap.ag <- df.heatmap.ag[,c(61,1:60)] # put row sums in front for convenience

low_paths <- df.heatmap.ag[df.heatmap.ag$sums < 4,]
df.heatmap.ag <- df.heatmap.ag[!df.heatmap.ag$sums < 4,]

df.heatmap.ag <- df.heatmap.ag[,c(-1)]

# Set annotation colors
ann_cols <- list(
  orig.ident = c(RCC1N = 'khaki1',RCC1T1 = 'skyblue',RCC1T2 = 'red')
)

# Plot heatmap
rcc1.heatmap <- pheatmap(as.matrix(df.heatmap.ag),
                        cluster_rows = F,
                        cluster_cols = F,
                        annotation_col = df.heatmap.anno,
                        annotation_colors = ann_cols,
                        show_colnames = F,
                        gaps_col = c(seq(3,60,by=3)),
                        main = 'RCC1 Pathways')
