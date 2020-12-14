setwd("E:/project/scripts")
library(Seurat)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)

#read data in and filter for oesophageal gland cells
seurat_object <- readRDS("data/GSE146736_adult_scseq_seurat.rds")
og <- seurat_object@active.ident
ogdf <- data.frame(og)
setDT(ogdf, keep.rownames = TRUE)[]
colnames(ogdf) <- c("cell", "tissue")
ogdf$tissue <- as.factor(ogdf$tissue)
d <- filter(ogdf, ogdf$tissue == "oesophageal gland")
c <- d$cell

#subset the main data object for the oesophageal cells
seurat_sub <- subset(seurat_object, cells = c)

#find top10 enriched genes for each slot in oesophageal subset
dat <- seurat_sub@assays[["RNA"]]@data
dat <- data.frame(dat)
datmean <- rowMeans(dat)
datmean <- data.frame(datmean)
setDT(datmean, keep.rownames = TRUE)[]
datmeansorted <- datmean[order(-datmean),]
head(datmeansorted, n=10)
# rn  datmean
# 1: Smp-040130 4.542484
# 2: Smp-072330 4.173694
# 3: Smp-335990 4.088275
# 4: Smp-155060 3.773790
# 5: Smp-021140 3.714004
# 6: Smp-082240 3.702959
# 7: Smp-031310 3.686800
# 8: Smp-032950 3.625991
# 9: Smp-030370 3.566331
# 10: Smp-098960 3.564980
dat <- seurat_sub@assays[["RNA"]]@counts
dat <- data.frame(dat)
datmean <- rowMeans(dat)
datmean <- data.frame(datmean)
setDT(datmean, keep.rownames = TRUE)[]
datmeansorted <- datmean[order(-datmean),]
head(datmeansorted, n=10)
# rn  datmean
# 1: Smp-040130 37.77966
# 2: Smp-072330 25.15254
# 3: Smp-335990 22.33898
# 4: Smp-032950 19.22034
# 5: Smp-082240 17.64407
# 6: Smp-095360 16.59322
# 7: Smp-155060 15.55932
# 8: Smp-066990 15.49153
# 9: Smp-021140 15.30508
# 10: Smp-335630 15.22034
dat <- seurat_sub@assays[["integrated"]]@data
dat <- data.frame(dat)
datmean <- rowMeans(dat)
datmean <- data.frame(datmean)
setDT(datmean, keep.rownames = TRUE)[]
datmeansorted <- datmean[order(-datmean),]
head(datmeansorted, n=10)
# rn  datmean
# 1: Smp-072330 4.215462
# 2: Smp-335990 4.019389
# 3: Smp-040130 3.939919
# 4: Smp-900060 3.925536
# 5: Smp-155060 3.730306
# 6: Smp-003770 3.549917
# 7: Smp-056970 3.530971
# 8: Smp-021140 3.525944
# 9: Smp-900090 3.497808
# 10: Smp-082240 3.425961

####THIS IS SLOT OF INTEREST ##### 
dat <- seurat_sub@assays[["integrated"]]@scale.data
dat <- data.frame(dat)
datmean <- rowMeans(dat)
datmean <- data.frame(datmean)
setDT(datmean, keep.rownames = TRUE)[]
datmeansorted <- datmean[order(-datmean),]
top10 <- head(datmeansorted, n=10)
# rn  datmean
# 1: Smp-172180 7.975959
# 2: Smp-124000 7.700045
# 3: Smp-328980 7.445943
# 4: Smp-085840 7.295113
# 5: Smp-125320 7.194271
# 6: Smp-010550 7.110686
# 7: Smp-123200 6.752361
# 8: Smp-176020 6.103522
# 9: Smp-315060 6.054716
# 10: Smp-171190 6.053953
FeaturePlot(object = seurat_object, cells = c, features = top10$rn, cols = c("orange", "blue"))

#Plots
library(ggplot2)
library(cowplot)
#map of whole dimplot
DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + theme(legend.position = "none") 
#highlight cells labelled as oesophageal gland amongst all cells
DimPlot(object = seurat_object, reduction = 'umap', cells.highlight = c) + theme(legend.position = "none") 
#plot only oesophageal gland cells 
DimPlot(object = seurat_object, cells = c, reduction = 'umap')

#recluster attempt
DimPlot(object = seurat_object, label = TRUE, pt.size = 0.5) + theme(legend.position = "none")

seurat_object[["ClusterNames_0.6"]] <- Idents(object = seurat_object)
seurat_object <- FindClusters(object = seurat_object, resolution = 1)
plot1 <- DimPlot(object = seurat_object, label = TRUE) + theme(legend.position = "none")
plot2 <- DimPlot(object = seurat_object, group.by = "ClusterNames_0.6", label = TRUE) + theme(legend.position = "none")
plot_grid(plot1, plot2)
plot3 <- DimPlot(object = seurat_object, label = TRUE, cells = c) 
plot4 <- DimPlot(object = seurat_object, group.by = "ClusterNames_0.6", label = TRUE, cells = c)
plot_grid(plot3, plot4)
DimPlot(object = seurat_object, reduction = 'umap', cells.highlight = c, label = TRUE) + theme(legend.position = "none") 
og.markers <- FindMarkers(object = seurat_object, ident.1 = 0, ident.2 = 54)
setDT(og.markers, keep.rownames = TRUE)[]

#get 0.05 as bonferoni corrected value 
dim(og.markers)
(0.05)/1293
#3.866976e-05
#filter by adjusted p value
filtered.marker <- filter(og.markers, p_val_adj <= 3.866976e-05)
ordered.marker <- filtered.marker[order(-avg_logFC)]
top5marker <- head(ordered.marker, n = 5)
ordered.marker <- filtered.marker[order(avg_logFC)]
bottom5marker <- head(ordered.marker, n = 5)
markers <- rbind(top5marker, bottom5marker)
#put results of markers into wormbase and saved as top10 foldchange
#plot for genes with largest fold changes 
FeaturePlot(object = seurat_object, cells = c, features = markers$rn, cols = c("orange", "blue"))

#Expression of three (MEGs 12, 16, 17) was confirmed in the anterior gland and 
#five (MEGs 8.1, 9, 11, 15 and 22) in the posterior gland.
#plot for genes shown to distinguish anterioir and posterior glands from literatire
MEGnames <- read.csv("ant_post_MEGs.csv")
colnames(MEGnames) <- c('name', 'ID')
num <- 1:8
for (i in num) {
  assign(paste('a', i, sep=''),FeaturePlot(object = seurat_object, cells = c, features = MEGnames$ID[i], cols = c("orange", "blue")) +
           labs(title = paste(MEGnames$name[i],"\n", MEGnames$ID[i])))
}
library(cowplot)
plot_grid(a1, a2, a3, a4, a5, a6, a7, a8)

#reluster attempt on og gland only ## desnt work
# DimPlot(object = seurat_sub, label = TRUE, pt.size = 0.5) + theme(legend.position = "none")
# seurat_sub[["ClusterNames_0.6"]] <- Idents(object = seurat_sub)
# seurat_sub <- FindClusters(object = seurat_sub, resolution = 0.8)
# plot1 <- DimPlot(object = seurat_sub, label = TRUE) + theme(legend.position = "none")
# plot2 <- DimPlot(object = seurat_sub, group.by = "ClusterNames_0.6", label = TRUE) + theme(legend.position = "none")
# plot_grid(plot1, plot2)

