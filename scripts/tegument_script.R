#Tegument analysis
setwd("E:/project/scripts")
library(Seurat)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
seurat_object <- readRDS("data/GSE146736_adult_scseq_seurat.rds")

teg <- seurat_object@active.ident
tegdf <- data.frame(teg)
setDT(tegdf, keep.rownames = TRUE)[]
colnames(tegdf) <- c("cell", "tissue")
tegdf$tissue <- as.factor(tegdf$tissue)
teg1 <- filter(tegdf, tegdf$tissue == "tegument 1")
teg2 <- filter(tegdf, tegdf$tissue == "tegument 2")
# teg3 <- filter(tegdf, tegdf$tissue == "early tsp-2+")
# teg4 <- filter(tegdf, tegdf$tissue == "egc+")
# teg5 <- filter(tegdf, tegdf$tissue == "meg-1+")
# teg6 <- filter(tegdf, tegdf$tissue == "zfp-1-1+")
# teg7 <- filter(tegdf, tegdf$tissue == "sm13+")

teg <- rbind(teg1, teg2)#, teg3, teg4, teg5, teg6, teg7)
c <- teg$cell


# all tegument progenitor cells ("early tsp-2+", "egc+", "meg-1+", "zfp-1-1+", and "sm13+") into a
# single cluster ("Tegument Progenitors"),

#subset the larger seurat object into a smaller one just for tegument
seurat_teg <- subset(seurat_object, cells = c)
#map of all ce;;s
DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + theme(legend.position = "none") 

#highlight cells labelled as tegument amongst all cells
DimPlot(object = seurat_object, reduction = 'umap', cells.highlight = c) + theme(legend.position = "none") 
#plot only tegument cells, this shows how the two tegument identities locate
DimPlot(object = seurat_object, cells = c, reduction = 'umap')


#extract coordinates of teg cells
DimPlot(object = seurat_teg, cells = c, reduction = 'umap')

coord <- seurat_teg[["umap"]]@cell.embeddings
coord2 <- as.data.frame(coord)
#make rownames into a column
library(data.table)
setDT(coord2, keep.rownames = TRUE)[]

#remove single cell outlier ACTATCTCACGTTGGC_9
#coord2 <- coord2[-41,]
#separate into upper and lower clusters, label, then combine
library(dplyr)
upper2 <- filter(coord2, (UMAP_2>=0))
upper2 <- filter(upper2, (UMAP_1<=0))
cluster <- c('upper2')
upper2 <- cbind(upper2, cluster)

lower1 <- filter(coord2, (UMAP_2<=0))
cluster <- c('lower1')
lower <- cbind(lower1, cluster)

upper1 <- filter(coord2, (UMAP_2>=0))
upper1 <- filter(upper1, (UMAP_1>=0))
cluster <- c('upper1')
upper1 <- cbind(upper1, cluster)

cluster_df <- rbind(lower, upper1, upper2)
colnames(cluster_df) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster')
cluster <- c('upper1', 'upper2', 'lower1')

# #add true neoblasts into a copy of the dataframe 
# neo <- seurat_object@active.ident
# neodf <- data.frame(neo)
# setDT(neodf, keep.rownames = TRUE)[]
# colnames(neodf) <- c("cell", "tissue")
# neodf$tissue <- as.factor(neodf$tissue)
# d <- filter(neodf, neodf$tissue == "neoblast")
# c <- d$cell
# seurat_teg_neo <- subset(seurat_object, cells = c)  
# coord <- seurat_teg_neo[["umap"]]@cell.embeddings
# coord2 <- as.data.frame(coord)
# library(data.table)
# setDT(coord2, keep.rownames = TRUE)[]
# cluster <- c('neoblast')
# neo_cells <- cbind(coord2, cluster)
# colnames(neo_cells) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster')
# neo_og_cells <- rbind(cluster_df, neo_cells)

# ### are the otegument cells not neoblasts
# library(tidyr)
# #tegument marker Sm25 Smp-346900
# #neoblast marker nanos2 = smp-051920
# dat2 <- seurat_teg@assays[["integrated"]]@scale.data
# dat3 <- seurat_teg_neo@assays[["integrated"]]@scale.data
# dat2 <- data.frame(dat2)
# dat3 <- data.frame(dat3)
# 
# newdat2 <- dat2['Smp-346900',]
# newdat3 <- dat3['Smp-346900',]
# 
# cn2 <- colnames(newdat2)
# cn3 <- colnames(newdat3)
# 
# setDT(newdat2, keep.rownames = TRUE)[]
# setDT(newdat3, keep.rownames = TRUE)[]
# 
# ogcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
# nbcountdf <- pivot_longer(newdat3, cols = cn3, names_to = "cell")
# 
# colnames(ogcountdf) <- c('ID', 'cell', 'scaled_counts_per_cell')
# colnames(nbcountdf) <- c('ID', 'cell', 'scaled_counts_per_cell')
# 
# ognbcountdf <- rbind(ogcountdf, nbcountdf)
# 
# nc <- ognbcountdf$cell
# 
# ogandneo <- merge(ognbcountdf, neo_og_cells, by = 'cell')
# 
# nc <- neo_og_cells$cell
# 
# 
# 
# ggplot(data = ogandneo, aes(x = cluster, y = scaled_counts_per_cell)) +
#   geom_violin() +
#   geom_jitter(width=0.15, alpha=0.5)
# 
# #means for each cluster group
# ogandneo %>% 
#   group_by(cluster) %>% 
#   summarise(mean(scaled_counts_per_cell))
# 
# 
# nbcountsummary <-  ogandneo %>%
#   group_by(cluster) %>%
#   summarise(mean = mean(scaled_counts_per_cell),
#             std = sd(scaled_counts_per_cell),
#             n = length(scaled_counts_per_cell),
#             se = std/sqrt(n))
# 
# #needs to be anova to compare the three clusters 
# mod <- aov(scaled_counts_per_cell ~ cluster, data = ogandneo)
# summary(mod)
# TukeyHSD(mod)
# plot(mod, 1)
# plot(mod, 2)
# #If data doesnt look normally distributed so use non-parametric test 
# 
# kw <- kruskal.test(scaled_counts_per_cell ~ cluster, data = ogandneo)
# pairwise.wilcox.test(ogandneo$scaled_counts_per_cell, ogandneo$cluster,
#                      p.adjust.method = "BH")




##### investigate surface proteins #######











## heatmap for all og proteins 
c <- cluster_df$cell 
meta <- cluster_df
names <- meta[,1]
library(tidyr)
meta$cell <- as.factor(meta$cell)
meta$cluster <- as.factor(meta$cluster)
meta <- data.frame(meta)
row.names(meta) <- meta$cell
cluster.info <- meta 
seurat_teg <- AddMetaData(object = seurat_teg, metadata = cluster.info)
DimPlot(object = seurat_teg, reduction = 'umap', group.by = "cluster")


library(RColorBrewer)
#by default heatmap looks as scale.data slot, add cells = c to remove outlier column
genes <- read.csv("gene list Master v2.csv")
colnames(genes) <- c('type','protein','gene_ID','table','tableName','tableID','genestableID','geneName','description','ignore')

genes2 <- filter(genes, (table==2))
genes3 <- filter(genes, (table==3))
genesteg <- rbind(genes2, genes3)

teg_short <- genesteg
rownames(teg_short) <- teg_short$gene_ID
teg_NA <- c('Smp-302170', 'Smp-121990', 'Smp-046290', 'Smp-194920', 'Smp-091240', 'Smp-142450', 'Smp-071250', 'Smp-082810', 'Smp-062300', 'Smp-305460', 'Smp-163750', 'Smp-074140', 'Smp-045500', 'Smp-017430', 'Smp-315900', 'Smp-311520', 'Smp-104500', 'Smp-091650', 'Smp-015020', 'Smp-344400', 'Smp-244150', 'Smp-048230', 'Smp-037540', 'Smp-022990', 'Smp-313560', 'Smp-137410', 'Smp-136690')
tegcleaned <- teg_short[!teg_short$gene_ID %in% teg_NA, ]

#scale_y_discrete labels from opposite way to plotting so label from an inverted dataframe
library(purrr)
tegflipped <- tegcleaned
tegflipped <- tegflipped %>% map_df(rev)


DoHeatmap(object = seurat_teg, features = tegcleaned$gene_ID, group.by = "cluster", draw.lines = TRUE, assay = 'RNA', slot = "counts")+ 
  scale_fill_gradientn(colours = c('black', 'orange')) +
  scale_y_discrete(labels = tegflipped$protein)


#try making heatmap read raw count expression slot

genestegflip <- genesteg
genestegflip <- genestegflip %>% map_df(rev)

DoHeatmap(object = seurat_teg, features = genesteg$gene_ID, group.by = "cluster", assay = 'RNA', slot = "counts") +
  scale_fill_gradientn(colours = c('black', 'orange')) +
  scale_y_discrete(labels = genestegflip$protein)




### RECLUSTER TEGUMENT 


DefaultAssay(object = seurat_teg) <- "RNA"
seurat_teg <- DietSeurat(seurat_teg,counts = TRUE, data = TRUE ,scale.data = FALSE,features = NULL, assays = 'RNA', dimreducs = NULL, graphs = NULL)

#follow Satija pre-clustering analysis 
seurat_teg <- NormalizeData(seurat_teg, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_teg <- FindVariableFeatures(seurat_teg, assay = 'RNA')
all.genes <- rownames(seurat_teg)
seurat_teg <- ScaleData(seurat_teg, features = all.genes)
#set npcs to be less than the number of cells (n = 46)
seurat_teg <- RunPCA(seurat_teg, npcs = 647,features = VariableFeatures(object = seurat_teg))
#check 'dimensionality' of data
seurat_teg <- JackStraw(seurat_teg, num.replicate = 100)
seurat_teg <- ScoreJackStraw(seurat_teg, dims = 1:20)
JackStrawPlot(seurat_teg, dims = 1:20)
#pc 1-20 look better than rest
ElbowPlot(seurat_teg)


#set dim = 1:25 because pc 1-20 looked really good so up 25 is probably good
seurat_teg <- FindNeighbors(seurat_teg, dims = 1:25)
#resolution 1 
seurat_teg <- FindClusters(seurat_teg, resolution = 1)


head(Idents(seurat_teg), 5)
#run UMAP
seurat_teg <- RunUMAP(seurat_teg, dims = 1:25)
#view
DimPlot(seurat_teg, reduction = "umap", pt.size = 3) 

#only run if you want to see how original clusters map to new clusters
# c <- cluster_df$cell 
# meta <- cluster_df
# names <- meta[,1]
# library(tidyr)
# meta$cell <- as.factor(meta$cell)
# meta$cluster <- as.factor(meta$cluster)
# meta <- data.frame(meta)
# row.names(meta) <- meta$cell
# cluster.info <- meta 
# seurat_teg <- AddMetaData(object = seurat_teg, metadata = cluster.info)
# 
# DimPlot(seurat_teg, reduction = "umap", pt.size = 3, group.by = "cluster" ) 

#find top 10 or 20 genes that drive the subclusters (i.e. top DE genes) and plots these as a heat map
top20 <- head(VariableFeatures(seurat_teg), 20)
plot1 <- VariableFeaturePlot(seurat_teg)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

top20_names <- read.csv('top20_tegument.csv')
colnames(top20_names) <- c('feature', 'name')
top20flipped <- top20_names %>% map_df(rev)


DoHeatmap(object = seurat_teg, features = top20, group.by = "seurat_clusters", assay = 'RNA', slot = "counts") +
  scale_fill_gradientn(colours = c('black', 'orange')) +
  scale_y_discrete(labels = top20flipped$name)

#heatmap for tegument proteins on reclustered cells
DoHeatmap(object = seurat_teg, features = genesteg$gene_ID, group.by = "seurat_clusters", assay = 'RNA', slot = "counts") +
  scale_fill_gradientn(colours = c('black', 'orange')) +
  scale_y_discrete(labels = genestegflip$protein)


## plot teg proteins on reculster UMAP
genesteg
num <- 1:66
for (i in num) {
  assign(paste('a', i, sep=''),FeaturePlot(object = seurat_teg, features = genesteg$gene_ID[i], cols = c("orange", "blue"), slot = "scale.data") +
           labs(title = paste(genesteg$protein[i],"\n", genesteg$gene_ID[i]))+
           theme(plot.title = element_text(size=8)))
}

teggraph <- ggarrange(a1,  a2,  a3,  a4,  a5,  a6,  a7,  a8,  a9,  a10, a11, a12, a13, a14,
                 a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28,
                 a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42,
                 a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56,
                 a57, a58, a59, a60, a61, a62, a63, a64, a65, a66, ncol=3, nrow=3)

png("tegument_UMAPs.png", width = 700, height =700)
plot(teggraph$`1`)
dev.off()
png("tegument_UMAPsb.png", width = 700, height =700)
plot(teggraph$`2`)
dev.off()
png("tegument_UMAPsc.png", width = 700, height =700)
plot(teggraph$`3`)
dev.off()
png("tegument_UMAPsd.png", width = 700, height =700)
plot(teggraph$`4`)
dev.off()
png("tegument_UMAPse.png", width = 700, height =700)
plot(teggraph$`5`)
dev.off()
png("tegument_UMAPsf.png", width = 700, height =700)
plot(teggraph$`6`)
dev.off()
png("tegument_UMAPsg.png", width = 700, height =700)
plot(teggraph$`7`)
dev.off()
png("tegument_UMAPsh.png", width = 700, height =700)
plot(teggraph$`8`)
dev.off()

