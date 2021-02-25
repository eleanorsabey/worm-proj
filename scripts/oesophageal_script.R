#Script for Oesophageal gland work; includes analysis of original clustering and reclustering.

#Load packages
setwd("E:/project/scripts")
library(Seurat)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RColorBrewer)
library(purrr)
library(writexl)
library(car)
library(ggsignif)
library(ggpubr)
library(grid)
library(ggplotify)


#Read Seurat object in and filter for oesophageal gland cells:
seurat_object <- readRDS("data/GSE146736_adult_scseq_seurat.rds")
og <- seurat_object@active.ident
ogdf <- data.frame(og)
setDT(ogdf, keep.rownames = TRUE)[]
colnames(ogdf) <- c("cell", "tissue")
ogdf$tissue <- as.factor(ogdf$tissue)
d <- filter(ogdf, ogdf$tissue == "oesophageal gland")
c <- d$cell

#Subset the main Seurat object for the oesophageal cells:
seurat_sub <- subset(seurat_object, cells = c)
all.genes <- rownames(seurat_sub)
seurat_sub <- ScaleData(seurat_sub, features = all.genes)


#Plot
#Map of whole Seurat object:
DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + theme(legend.position = "none") 
#Oesophageal gland cells highlighted amongst all cells:
DimPlot(object = seurat_object, reduction = 'umap', cells.highlight = c) + theme(legend.position = "none") 
#Plot only oesophageal gland cells:
DimPlot(object = seurat_object, cells = c, reduction = 'umap')

#Investigation into oesophageal cells; includes comparison of sub clusters with neoblast cells.

#Extract coordinates of oesophageal cells:
DimPlot(object = seurat_sub, cells = c, reduction = 'umap')

coord <- seurat_sub[["umap"]]@cell.embeddings
coord2 <- as.data.frame(coord)
setDT(coord2, keep.rownames = TRUE)[]
#remove single cell outlier ACTATCTCACGTTGGC_9
coord2 <- coord2[-41,]
#separate into upper and lower clusters, label, then combine
upper <- filter(coord2, (UMAP_2>=0))
cluster <- c('upper')
upper <- cbind(upper, cluster)
lower <- filter(coord2, (UMAP_2<=0))
cluster <- c('lower')
lower <- cbind(lower, cluster)
cluster_df <- rbind(lower, upper)
colnames(cluster_df) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster')

#Add true neoblasts into a copy of the dataframe: 
neo <- seurat_object@active.ident
neodf <- data.frame(neo)
setDT(neodf, keep.rownames = TRUE)[]
colnames(neodf) <- c("cell", "tissue")
neodf$tissue <- as.factor(neodf$tissue)
d <- filter(neodf, neodf$tissue == "neoblast")
c <- d$cell
seurat_sub_neo <- subset(seurat_object, cells = c) 
coord <- seurat_sub_neo[["umap"]]@cell.embeddings
coord2 <- as.data.frame(coord)
setDT(coord2, keep.rownames = TRUE)[]
cluster <- c('neoblast')
neo_cells <- cbind(coord2, cluster)
colnames(neo_cells) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster')
neo_og_cells <- rbind(cluster_df, neo_cells)

#make a seurat object with Og and Neo cells, then re scale.
allcells <- neo_og_cells$cell
seurat_og_neo <- subset(seurat_object, cells = allcells)
all.genes <- rownames(seurat_og_neo)
seurat_og_neo <- ScaleData(seurat_og_neo, features = all.genes)

#Compare expression of any gene between the sub-clusters of oesophageal cells and the neoblast cluster:
#Genes used in analysis: 
#MEG-8.2 = Smp-172180 this ID was used to highlight OG cells in the Wendt Paper.
#MEG-4.2 = Smp-085840 and MEG-14 = Smp-124000, meant to mark OG gland (https://doi.org/10.1371/journal.pntd.0002337).
#Neoblast marker nanos2 = Smp-051920 this ID was used to highlight OG cells in the Wendt Paper.
dat2 <- seurat_og_neo@assays[["integrated"]]@scale.data
dat2 <- data.frame(dat2)


#these genes dont seem to be working
genes <- c('Smp-051920','Smp-105360','Smp-175590','Smp-125320','Smp-010550', 'Smp-172180', 'Smp-085840', 'Smp-124000')
names <- c('nanos-2','notch','fgfra','meg-9','meg-15', 'meg-8.2', 'meg-4.2','meg-14')
stat <- cbind(genes, names)
stat <- data.frame(stat)
num <- 1:8

fileConn<-file("output.txt")
sink(fileConn, append=TRUE)

for (i in num){
  
  n <- stat$names[i]
  id <- stat$genes[i]
  
  print(paste('Stats for gene', n, id))
 
  newdat2 <- dat2[id,]

  cn2 <- colnames(newdat2)

  setDT(newdat2, keep.rownames = TRUE)[]

  ognbcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")

  colnames(ognbcountdf) <- c('ID', 'cell', 'expression')

  ogandneo <- merge(ognbcountdf, neo_og_cells, by = 'cell')
  
  #Mean expression for each group:
  print('Mean expression:')
  print(ogandneo %>% 
          group_by(cluster) %>% 
          summarise(mean(expression)))  
  
  my_comparisons = list(c("lower","neoblast"),c("upper", "lower"), c("upper", "neoblast"))
  assign(paste("Plot", i, sep = ''),ggplot(data = ogandneo, aes(x = cluster, y = expression)) +
    geom_violin() +
    labs(title = paste(n,"\n", id))+
      theme(plot.title = element_text(size=8))+
    geom_jitter(width=0.15, alpha=0.5) +
    stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = 'p.signif', paired = FALSE))
  
  #Anova to compare expression across the three clusters:
  mod <- aov(expression ~ cluster, data = ogandneo)
  #View anova:
  print('Parametric ANOVA and TukeyHSD results:')
  print(summary(mod))
  print(TukeyHSD(mod))
  
  #Check normality:
  assign(paste("hist", i, sep = ''), as.ggplot(~hist(mod$residuals)))
  #Check equal varience:
  assign(paste("box", i, sep = ''), as.ggplot(~boxplot(expression ~ cluster, data = ogandneo)))
  

  print("Non-parametric - Kruskal Wallis and unpaired Wilcox test:")
  #If varience isnt equal: Welch's anova 
  print(kruskal.test(expression ~ cluster, data = ogandneo))
  print(pairwise.wilcox.test(ogandneo$expression, ogandneo$cluster,
                       p.adjust.method = "none", paired = FALSE))
  
}
close(fileConn)
sink()
statplot <- ggarrange(Plot1, hist1, box1, Plot2, hist2, box2, Plot3, hist3, box3, Plot4,  hist4, box4, ncol=3, nrow=4)
statplot2 <- ggarrange(Plot5, hist5, box5, Plot6,  hist6, box6, Plot7, hist7, box7, Plot8, hist8, box8, ncol=3, nrow=4)
plot(statplot)
png("statsplotOG1.png", width = 470, height = 800)
plot(statplot)
dev.off()
png("statsplotOG2.png", width = 470, height = 800)
plot(statplot2)
dev.off()
  
  




#Further analysis by investigating expression of associated oesophageal proteins across the sub-clusters of cells:
#Add sub-cluster labelling to Seurat object:
c <- cluster_df$cell 
meta <- cluster_df[,-c(2:3)]
cell <- "ACTATCTCACGTTGGC_9"
cluster <- "singleton"
singleton <- data.frame(cell, cluster)
meta<- rbind(meta, singleton)
names <- meta[,1]
meta$cell <- as.factor(meta$cell)
meta$cluster <- as.factor(meta$cluster)
meta <- data.frame(meta)
row.names(meta) <- meta$cell
cluster.info <- meta 
seurat_sub <- AddMetaData(object = seurat_sub, metadata = cluster.info)
#Check if labelling worked:
DimPlot(object = seurat_sub, reduction = 'umap', group.by = "cluster")

#Load protein data:
eso_proteins <- read.csv("esophageal_gland_proteins.csv")
colnames(eso_proteins) <- c('name', 'ID')


#Scale_y_discrete labels the plot in reverse, so we need to make an inverted dataframe:
eso_proteinsflip <- eso_proteins
eso_proteinsflip <- eso_proteinsflip %>% map_df(rev)


#Heatmap using scale.data
DoHeatmap(object = seurat_sub, group.by= 'cluster', features = eso_proteins$ID, draw.lines = TRUE)+ 
    scale_fill_gradientn(colours = c('blue', 'black', 'orange')) +
  scale_y_discrete(labels = eso_proteinsflip$name)
  

#It is decidede that the upper cluster resembles oesophageal cells more than neoblast so will not be removed from the data before reclustering. 
#The outlier cell is also not removed becase heatmap expression fits with lower cluster. 

#Re-cluster the lower cluster and assess differencial expression that may indicate anterior/posterior tissues.



#Set default assay to RNA and remove the integrated assay:
DefaultAssay(object = seurat_sub) <- "RNA"
seurat_sub <- DietSeurat(seurat_sub,counts = TRUE, data = TRUE ,scale.data = FALSE,features = NULL, assays = 'RNA', dimreducs = NULL, graphs = NULL)

#Re-cluster following Satija pre-clustering analysis tutorial:
seurat_sub <- NormalizeData(seurat_sub, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_sub <- FindVariableFeatures(seurat_sub, assay = 'RNA')
all.genes <- rownames(seurat_sub)
seurat_sub <- ScaleData(seurat_sub, features = all.genes)
#Set npcs to be less than the number of cells (n = 59)
seurat_sub <- RunPCA(seurat_sub, npcs = 58,features = VariableFeatures(object = seurat_sub))
#Check 'dimensionality' of data, to know how mnay PCs to use in clustering functions:
seurat_sub <- JackStraw(seurat_sub, num.replicate = 100)
seurat_sub <- ScoreJackStraw(seurat_sub, dims = 1:20)
JackStrawPlot(seurat_sub, dims = 1:20)
ElbowPlot(seurat_sub)
#Quality drops after PC 1-4, John advised to go anywhere up to 10.
#Set dim = 1:4
seurat_sub <- FindNeighbors(seurat_sub, dims = 1:4)
seurat_sub <- FindClusters(seurat_sub, resolution = 1)
head(Idents(seurat_sub), 5)
seurat_sub <- RunUMAP(seurat_sub, dims = 1:4)
#View new clustering:
DimPlot(seurat_sub, reduction = "umap", pt.size = 3, cols = c('#7AD169', '#928DE0', '#FAA637')) 

#ONly run to see how old clustering maps to new:
cluster.info <- meta
seurat_sub <- AddMetaData(object = seurat_sub, metadata = cluster.info)

DimPlot(seurat_sub, reduction = "umap", pt.size = 3, group.by = "cluster" )

#Heatmap to assess expression of oesophageal proteins on new clustering:
#This uses scale.data expression:
eso_proteinsflip <- eso_proteins
eso_proteinsflip <- eso_proteinsflip %>% map_df(rev)

DoHeatmap(object = seurat_sub, features = eso_proteins$ID, group.by = "seurat_clusters", assay = 'RNA', slot = "scale.data", group.colors = c('#7AD169', '#928DE0', '#FAA637')) +
  scale_fill_gradientn(colours = c('blue', 'black', 'orange'), breaks=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3),limits=c(-2.5,2.5)) +
  scale_y_discrete(labels = eso_proteinsflip$name) 

DoHeatmap(object = seurat_sub, features = eso_proteins$ID, group.by = "seurat_clusters", assay = 'RNA', slot = "data", group.colors = c('#7AD169', '#928DE0', '#FAA637')) +
  scale_fill_gradientn(colours = c('black', 'orange')) +
  scale_y_discrete(labels = eso_proteinsflip$name) 

#Compare expression of specific proteins across new clusters:

og <- seurat_lower@active.ident
ogdf <- data.frame(og)
setDT(ogdf, keep.rownames = TRUE)[]
colnames(ogdf) <- c("cell", "cluster")
ogdf$cluster <- as.factor(ogdf$cluster)
cluster_df <- ogdf

#Attempts to find markers of each cluster

# #1. VariableFeature approach:
# top20 <- head(VariableFeatures(seurat_sub), 20)
# plot1 <- VariableFeaturePlot(seurat_sub)
# plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
# #Enter top20 names into WormBase BioMart and find names.
# top20_names <- read.csv('top20OG.csv')
# colnames(top20_names) <- c('feature', 'name')
# top20flipped <- top20_names %>% map_df(rev)
# 
# #Plot top20 on Heatmap:
# DoHeatmap(object = seurat_sub, features = top20, group.by = "seurat_clusters", assay = 'RNA', slot = "scale.data", group.colors = c('#7AD169', '#928DE0', '#FAA637')) +
#   scale_fill_gradientn(colours = c('blue', 'black', 'orange')) +
#   scale_y_discrete(labels = top20flipped$name)


#2. Find markers that dinstunguish the sub clusters from each other:
#In general, an AUC of 0.5 suggests no discrimination (i.e., ability to diagnose patients with and without the disease or condition based on the test), 0.7 to 0.8 is considered acceptable, 0.8 to 0.9 is considered excellent, and more than 0.9 is considered outstanding.

allmarkers <- FindAllMarkers(seurat_sub, assay = 'RNA', test.use =  'roc', only.pos = TRUE, return.thresh = 0)
markersOG <- filter(allmarkers, myAUC >= 0.8)
write_xlsx(markersOG,"markersOG2.xlsx")

markersannotatedOG <- read.csv("markersOGannotated2.csv")

colnames(markersannotatedOG) <- c('myAUC',"avg_diff","power","avg_logFC","pct.1","pct.2","cluster","gene","name")
markersflipped <- markersannotatedOG %>% map_df(rev)

DoHeatmap(object = seurat_sub, features = markersOG$gene, group.by = "seurat_clusters", assay = 'RNA', slot = "scale.data", group.colors = c('#7AD169', '#928DE0', '#FAA637')) +
  scale_fill_gradientn(colours = c('blue', 'black', 'orange')) +
  scale_y_discrete(labels = markersflipped$name)

#3. Find markers that distinguish the sub clusters from all other clusters.
##Doesn't work on my laptop due to insufficient memory
test <- seurat_lower@active.ident
testdf <- data.frame(test)
setDT(testdf, keep.rownames = TRUE)[]
colnames(testdf) <- c("cell", "cluster")
testdf$cluster <- as.factor(testdf$cluster)
d0 <- filter(testdf, testdf$cluster == "0")
c0 <- d0$cell
d1 <- filter(testdf, testdf$cluster == "1")
c1 <- d1$cell
d2 <- filter(testdf, testdf$cluster == "2")
c2 <- d2$cell

testmarkers <- FindMarkers(seurat_object, ident.1 = c0, ident.2 = NULL)
testmarkers <- FindMarkers(seurat_object, ident.1 = c1, ident.2 = NULL)
testmarkers <- FindMarkers(seurat_object, ident.1 = c2, ident.2 = NULL)






