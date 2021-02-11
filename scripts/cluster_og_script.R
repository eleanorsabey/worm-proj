setwd("E:/project/scripts")
library(Seurat)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggpubr)

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


coord <- seurat_sub[["umap"]]@cell.embeddings
coord2 <- as.data.frame(coord)
#make rownames into a column
library(data.table)
setDT(coord2, keep.rownames = TRUE)[]

#remove single cell outlier ACTATCTCACGTTGGC_9
coord2 <- coord2[-41,]
#separate into upper and lower clusters, label, then combine
library(dplyr)
upper <- filter(coord2, (UMAP_2>=0))
cluster <- c('upper')
upper <- cbind(upper, cluster)
lower <- filter(coord2, (UMAP_2<=0))
cluster <- c('lower')
lower <- cbind(lower, cluster)
cluster_df <- rbind(lower, upper)
colnames(cluster_df) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster')

colnames(lower) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster') 
lnames <- lower$cell
namelist <- as.list(lnames)
lapply(namelist, write, "lowercells.txt", append=TRUE)
my_data <- read.delim("lowercells.txt", header = F)



seurat_lower <- subset(seurat_sub, cells = lnames)

#view lower cluster
DimPlot(object = seurat_lower, label = FALSE, pt.size = 3, cols = 'gray') #+ theme(legend.position = "none")

#set default assay to RNA and remove integrated assay
DefaultAssay(object = seurat_lower) <- "RNA"
seurat_lower <- DietSeurat(seurat_lower,counts = TRUE, data = TRUE ,scale.data = FALSE,features = NULL, assays = 'RNA', dimreducs = NULL, graphs = NULL)

#follow Satija pre-clustering analysis 
seurat_lower <- NormalizeData(seurat_lower, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_lower <- FindVariableFeatures(seurat_lower, assay = 'RNA')
all.genes <- rownames(seurat_lower)
seurat_lower <- ScaleData(seurat_lower, features = all.genes)
#set npcs to be less than the number of cells (n = 46)
seurat_lower <- RunPCA(seurat_lower, npcs = 45,features = VariableFeatures(object = seurat_lower))
#check 'dimensionality' of data
seurat_lower <- JackStraw(seurat_lower, num.replicate = 200)
seurat_lower <- ScoreJackStraw(seurat_lower, dims = 1:20)
JackStrawPlot(seurat_lower, dims = 1:20)
#pc 1-4 look better than rest
ElbowPlot(seurat_lower)


#set dim = 1:4 because pc 1-4 look okay
seurat_lower <- FindNeighbors(seurat_lower, dims = 1:6)
#resolution 0.9 makes each cluster on UMAP it's own colour roughly
seurat_lower <- FindClusters(seurat_lower, resolution = 1)

                            
head(Idents(seurat_lower), 5)
#run UMAP
seurat_lower <- RunUMAP(seurat_lower, dims = 1:6)
#view
DimPlot(seurat_lower, reduction = "umap", pt.size = 3, cols = c('#7AD169', '#928DE0', '#FAA637')) 

eso_proteins <- read.csv("esophageal_gland_proteins.csv")
colnames(eso_proteins) <- c('name', 'ID')
num <- 1:71
#work out what slot is plotted
#Expression of three (MEGs 12 Smp-152630, 16 Smp-158890, 17 Smp-180620) was confirmed in the anterior gland and 

#five (MEGs 8.1, 9, 11, 15 and 22) in the posterior gland.

for (i in num) {
  assign(paste('a', i, sep=''),FeaturePlot(object = seurat_lower, features = eso_proteins$ID[i], cols = c("orange", "blue"), slot = "data", pt.size = 3) +
           labs(title = paste(eso_proteins$name[i],"\n", eso_proteins$ID[i]))+
           theme(plot.title = element_text(size=8)))
}

plot(a1)

eso <- ggarrange(a1,  a2,  a3,  a4,  a5,  a6,  a7,  a8,  a9,  a10, a11, a12, a13, a14,
                 a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28,
                 a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42,
                 a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56,
                 a57, a58, a59, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a70,
                 a71, ncol=3, nrow=3)

png("esophageal_UMAPs.png", width = 805, height =700)
plot(eso$`1`)
dev.off()
png("esophageal_UMAPsb.png", width = 805, height =700)
plot(eso$`2`)
dev.off()
png("esophageal_UMAPsc.png", width = 805, height =700)
plot(eso$`3`)
dev.off()
png("esophageal_UMAPsd.png", width = 805, height =700)
plot(eso$`4`)
dev.off()
png("esophageal_UMAPse.png", width = 805, height =700)
plot(eso$`5`)
dev.off()
png("esophageal_UMAPsf.png", width = 805, height =700)
plot(eso$`6`)
dev.off()
png("esophageal_UMAPsg.png", width = 805, height =700)
plot(eso$`7`)
dev.off()
png("esophageal_UMAPsh.png", width = 805, height =700)
plot(eso$`8`)
dev.off()


cluster0.markers <- FindMarkers(seurat_lower, ident.1 = 0,ident.2 = c(1, 2), min.pct = 0.50, only.pos = T)
head(cluster0.markers, n = 5)

cluster1.markers <- FindMarkers(seurat_lower, ident.1 = 1, ident.2 = c(0, 2), min.pct = 0.50, only.pos = T)
head(cluster1.markers, n = 5)

cluster2.markers <- FindMarkers(seurat_lower, ident.1 = 2,ident.2 = c(1, 2), min.pct = 0.50, only.pos = T)
head(cluster2.markers, n = 5)
