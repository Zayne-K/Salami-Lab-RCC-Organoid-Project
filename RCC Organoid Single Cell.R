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
# library(metap)
# library(multtest)

setdir <- "/home/kzayne/Salami-Lab-RCC-Organoid-Project/"
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

### Cell markers file
kid.mrkrs <- read.csv("Final kidney markers.csv",header = T)
kid.mrkrs <- kid.mrkrs[!kid.mrkrs$cell_name %in% c("Neutrophil","Cancer stem cell"),]
mrkr.list <- as.list(as.character(kid.mrkrs$Symbol))
names(mrkr.list) <- kid.mrkrs$cell_name

for(i in 1:length(mrkr.list)){
  mrkr.list[[i]] <- unlist(strsplit(mrkr.list[[i]], "," ))
}

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
clusters <- DimPlot(rcc1, reduction = "umap",pt.size = 0.5,group.by = "seurat_clusters",label = T,repel = T,label.size = 5,order=T) + ggtitle("By UMAP Clusters")
clusters

ggarrange(ggarrange(condition,clusters,nrow=1),orig.stim,ncol=1)

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
rcc1$ca9 <- rcc1@assays$RNA$scale.data['CA9',]
rcc1$cellassign <- ifelse(rcc1$cellassign == 'Proximal tubular cell' & rcc1$ca9 > 0,'Proximal Tubular cell + CA9',rcc1$cellassign)

custom1 <- DimPlot(rcc1, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom1

custom2 <- DimPlot(rcc1, reduction = "umap", label = TRUE, repel = TRUE,
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom2

cairo_pdf(paste0(setdir,"RCC1_post_harmony.pdf"),width = 22,height = 20)
print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
dev.off()

cairo_pdf(paste0(setdir,"RCC1 UMAP cell clusters.pdf"),width = 22,height = 11)
print(custom1)
dev.off()

FeaturePlot(rcc1,
            repel=F,order = T,
            features = c("VEGFA","CA9","NCAM1","EPCAM"),
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

custom2 <- DimPlot(rcc2, reduction = "umap", label = TRUE, repel = TRUE,
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom2

cairo_pdf(paste0(setdir,"RCC2_post_harmony.pdf"),width = 22,height = 20)
print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
dev.off()

cairo_pdf(paste0("RCC2 UMAP cell clusters.pdf"),width = 22,height = 11)
print(custom1)
dev.off()

FeaturePlot(rcc2,
            repel=F,order = T,
            features = c("VEGFA","CA9","NCAM1","EPCAM"),
            split.by = "orig.ident",pt.size = 0.5)
#####

# Combine all seurat objects for RCC3
#####
set.seed(555)
rcc3 <- merge(rcc3n, y = c(rcc3t1,rcc3t2), add.cell.ids = c("RCC3N", "RCC3T1", "RCC3T2")) %>%
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

custom1 <- DimPlot(rcc3, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom1

custom2 <- DimPlot(rcc3, reduction = "umap", label = TRUE, repel = TRUE,
                   group.by = 'cellassign',pt.size = 0.5,label.size = 5,order = T) +
  ggtitle("By cell type")
custom2

cairo_pdf(paste0(setdir,"RCC3_post_harmony.pdf"),width = 22,height = 20)
print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
dev.off()

cairo_pdf(paste0(setdir,"RCC3 UMAP cell clusters.pdf"),width = 22,height = 11)
print(custom1)
dev.off()

FeaturePlot(rcc3,
            repel=F,order = T,
            features = c("VEGFA","CA9","NCAM1","EPCAM"),
            split.by = "orig.ident",pt.size = 0.5)

#####

# Combine all seurat objects for RCC1 RCC2 RCC3
#####
set.seed(555)
rcc.all <- merge(rcc1n, y = c(rcc1t1,rcc1t2,rcc2n,rcc2t1,rcc2t2,rcc3n,rcc3t1,rcc3t2),
                 add.cell.ids = c("RCC1N","RCC1T1","RCC1T2","RCC2N","RCC2T1","RCC2T2","RCC3N","RCC3T1","RCC3T2")) %>%
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
                              rcc.filt.split[[8]],rcc.filt.split[[9]]),
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

cairo_pdf(paste0(setdir,"rcc_all_harmony_plots.pdf"),width = 18,height = 13)
print(ggarrange(p1o,p1n,p2o,p2n,nrow=2,ncol=2))
dev.off()

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

cairo_pdf(paste0(setdir,"rcc_all_seurat_clusters.pdf"),width = 22,height = 20)
print(orig.stim)
dev.off()

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

cairo_pdf(paste0(setdir,"rcc_all_clusters.pdf"),width = 14,height = 11)
print(ggarrange(condition,clusters,nrow=1))
dev.off()

# Cell type assignment
# Assign clusters
# get cell-type by cell matrix
DefaultAssay(rcc.filt) <- "RNA"
rcc.filt <- ScaleData(rcc.filt,verbose = T,features = rownames(rcc.filt))
es.max = sctype_score(scRNAseqData = rcc.filt@assays$RNA$scale.data, 
                      scaled = TRUE, 
                      gs = mrkr.list)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(rcc.filt@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rcc.filt@meta.data[rcc.filt@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rcc.filt@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
rcc.filt$ca9 <- rcc.filt@assays$RNA$scale.data['CA9',]

### UMAP
rcc.filt@meta.data$cellassign = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  rcc.filt@meta.data$cellassign[rcc.filt@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
rcc.filt$cellassign <- ifelse(rcc.filt$cellassign == 'Proximal tubular cell' & rcc.filt$ca9 > 0,'Proximal Tubular cell + CA9',rcc.filt$cellassign)

custom1 <- DimPlot(rcc.filt, reduction = "umap", label = TRUE, repel = TRUE,
                   split.by = "orig.ident",ncol=3,
                   group.by = 'cellassign',pt.size = 0.3,label.size = 5,order = T) + 
  ggtitle("By cell type")
custom1

custom2 <- DimPlot(rcc.filt, reduction = "umap", label = TRUE, repel = TRUE,
                   group.by = 'cellassign',pt.size = 0.3,label.size = 5,order = T) + 
  ggtitle("By cell type")
custom2

cowplot::plot_grid(ncol = 2, DimPlot(rcc.filt, group.by = "orig.ident") + NoAxes(),
                   DimPlot(rcc.filt, group.by = "cellassign") + NoAxes())

cairo_pdf(paste0(setdir,"RCC all plots.pdf"),width = 22,height = 20)
print(ggarrange(ggarrange(condition,clusters,custom2,nrow=1),orig.stim,custom1,ncol=1))
dev.off()

cairo_pdf(paste0(setdir,"RCC all UMAP cell clusters.pdf"),width = 22,height = 15)
print(custom1)
dev.off()

cairo_pdf(paste0(setdir,"RCC all VEGFA CA9.pdf"),width = 30,height = 14)
print(FeaturePlot(rcc.filt,
                  repel = T,order = T,
                  features = c("VEGFA","CA9","EGFR"),
                  split.by = "orig.ident",
                  pt.size = 0.3))
dev.off()

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

cairo_pdf(paste0(setdir,"RCC all UMAP cell clusters by tissue type.pdf"),width = 18,height = 9)
print(custom3)
dev.off()

FeaturePlot(rcc.filt,
            repel = T,order = T,
            features = c("VEGFA","CA9","EGFR"),
            split.by = "tissue.type",
            cols = c("gray","red"),
            pt.size = 0.5)
#####

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
fc <- 1

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
hm.sym <- GSA.read.gmt("/home/kzayne/h.all.v7.0.symbols.gmt")    # download all available genesets genesymbols .gmt file from MSigDB
names(hm.sym$genesets) <- hm.sym$geneset.names
hm.sym <- hm.sym$genesets

# Xcell pathways
xcell.sig <- read.csv("/home/kzayne/XCell_pathways.csv",header = T)
xcell <- as.list(as.character(xcell.sig$Genes))
names(xcell) <- xcell.sig$Celltype_Source_ID

for(i in 1:length(xcell)){
  xcell[[i]] <- unlist(strsplit(xcell[[i]], "," ))
}

# RCC1
#####
# Identify DEGs between T1 and Normal
rcc1.T1vN <- FindMarkers(rcc.norm, ident.1 = "RCC1T1", ident.2 = "RCC1N")

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

#### ssgsea
library(scGSVA)
rcc1.norm <- subset(rcc.norm, subset = stim.new == 'RCC1')
rcc2.norm <- subset(rcc.norm, subset = stim.new == 'RCC2')
rcc3.norm <- subset(rcc.norm, subset = stim.new == 'RCC3')

rcc1.anno <- buildAnnot(species = 'human',keytype='SYMBOL',anntype='GO')
res <- scgsva(rcc3.norm, rcc1.anno,method = 'ssgsea')


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