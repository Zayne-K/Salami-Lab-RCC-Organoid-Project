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
