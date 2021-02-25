#Tegument analysis
setwd("E:/project/scripts")
library(Seurat)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(writexl)
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

#Add true neoblasts into a copy of the dataframe: 
neo <- seurat_object@active.ident
neodf <- data.frame(neo)
setDT(neodf, keep.rownames = TRUE)[]
colnames(neodf) <- c("cell", "tissue")
neodf$tissue <- as.factor(neodf$tissue)
d <- filter(neodf, neodf$tissue == "neoblast")
c <- d$cell
seurat_teg_neo <- subset(seurat_object, cells = c) 

coord <- seurat_teg_neo[["umap"]]@cell.embeddings
coord2 <- as.data.frame(coord)
setDT(coord2, keep.rownames = TRUE)[]
cluster <- c('neoblast')
neo_cells <- cbind(coord2, cluster)
colnames(neo_cells) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster')
neo_teg_cells <- rbind(cluster_df, neo_cells)

#make a seurat object with Og and Neo cells, then re scale.
allcells <- neo_teg_cells$cell
seurat_teg_neo <- subset(seurat_object, cells = allcells)
all.genes <- rownames(seurat_teg_neo)
seurat_teg_neo <- ScaleData(seurat_teg_neo, features = all.genes)

###############
dat2 <- seurat_teg_neo@assays[["integrated"]]@scale.data
dat2 <- data.frame(dat2)

#genes to be interested in:
genes <- c('Smp-051920','Smp-105360','Smp-175590','Smp-045200','Smp-346900')
names <- c('nanos-2','notch','fgfra','tal','sm25')
stat <- cbind(genes, names)
stat <- data.frame(stat)
num <- 1:8

fileConn<-file("outputTegstats.txt")
sink(fileConn, append=TRUE)

for (i in num){
  
  n <- stat$names[i]
  id <- stat$genes[i]
  
  print(paste('Stats for gene', n, id))
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  tegnbcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(tegnbcountdf) <- c('ID', 'cell', 'expression')
  
  tegandneo <- merge(tegnbcountdf, neo_teg_cells, by = 'cell')
  
  #Mean expression for each group:
  print('Mean expression:')
  print(tegandneo %>% 
          group_by(cluster) %>% 
          summarise(mean(expression)))  
  
  my_comparisons = list(c("lower1","neoblast"),c("upper1", "lower1"), c("upper1", "neoblast"), c("upper2", "neoblast"), c("upper2","upper1"), c("upper2","lower1"))
  assign(paste("Plot", i, sep = ''),ggplot(data = tegandneo, aes(x = cluster, y = expression)) +
           geom_violin() +
           labs(title = paste(n,"\n", id))+
           theme(plot.title = element_text(size=8))+
           geom_jitter(width=0.15, alpha=0.5) +
           stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = 'p.signif', paired = FALSE))
  
  #Anova to compare expression across the three clusters:
  mod <- aov(expression ~ cluster, data = tegandneo)
  #View anova:
  print('Parametric ANOVA and TukeyHSD results:')
  print(summary(mod))
  print(TukeyHSD(mod))
  
  #Check normality:
  assign(paste("hist", i, sep = ''), as.ggplot(~hist(mod$residuals)))
  #Check equal varience:
  assign(paste("box", i, sep = ''), as.ggplot(~boxplot(expression ~ cluster, data = tegandneo)))
  
  
  print("Non-parametric - Kruskal Wallis and unpaired Wilcox test:")
  #If varience isnt equal: Welch's anova 
  print(kruskal.test(expression ~ cluster, data = tegandneo))
  print(pairwise.wilcox.test(tegandneo$expression, tegandneo$cluster,
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
plot(Plot1)
###############
dat3 <- seurat_teg_neo@assays[["integrated"]]@scale.data
dat3 <- data.frame(dat3)

#Change Smp-ID for genes of interest 
newdat3 <- dat3['Smp-105360',]

cn3 <- colnames(newdat3)

setDT(newdat3, keep.rownames = TRUE)[]

tegcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
nbcountdf <- pivot_longer(newdat3, cols = cn3, names_to = "cell")

colnames(tegcountdf) <- c('ID', 'cell', 'expression')
colnames(nbcountdf) <- c('ID', 'cell', 'expression')

tegnbcountdf <- rbind(tegcountdf, nbcountdf)

nc <- tegnbcountdf$cell

tegandneo <- merge(tegnbcountdf, neo_teg_cells, by = 'cell')

nc <- neo_teg_cells$cell

#Anova to compare expression across the three clusters:
mod <- aov(expression ~ cluster, data = tegandneo)
#View anova:
summary(mod)
TukeyHSD(mod)

#Check normality:
hist(mod$residuals)
qqPlot(mod$residuals)
#Check equal varience:
boxplot(expression ~ cluster,
        data = tegandneo)

#If varience isnt equal: Welch's anova 
oneway.test(expression ~ cluster,
            data = tegandneo,
            var.equal = FALSE )

pairwise.wilcox.test(tegandneo$expression, tegandneo$cluster,
                     p.adjust.method = "none", paired = FALSE)


#If data doesnt look normally distributed so use non-parametric test 

kw <- kruskal.test(expression ~ cluster, data = ogandneo)


#Plot of expression
my_comparisons = list(c("lower1","neoblast"),c("upper1", "lower1"), c("upper2", "neoblast"), c('upper1', 'upper2'), c('upper1','neoblast' ))
ggplot(data = tegandneo, aes(x = cluster, y = expression)) +
  geom_violin() +
  geom_jitter(width=0.15, alpha=0.5) +
  stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = 'p.signif', paired = FALSE)

#Mean expression for each group:
ogandneo %>% 
  group_by(cluster) %>% 
  summarise(mean(expression))





## heatmap for all teg proteins 
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


DoHeatmap(object = seurat_teg, features = tegcleaned$gene_ID, group.by = "cluster", draw.lines = TRUE, assay = 'integrated', slot = "scale.data")+ 
  scale_fill_gradientn(colours = c('blue','black', 'orange')) +
  scale_y_discrete(labels = tegflipped$protein)


DoHeatmap(object = seurat_teg_neo, features = tegcleaned$gene_ID, draw.lines = TRUE, assay = 'integrated', slot = "scale.data")+ 
  scale_fill_gradientn(colours = c('blue','black', 'orange')) +
  scale_y_discrete(labels = tegflipped$protein)

#try making heatmap read raw count expression slot

genestegflip <- genesteg
genestegflip <- genestegflip %>% map_df(rev)

DoHeatmap(object = seurat_teg, features = genesteg$gene_ID, group.by = "cluster", assay = 'RNA', slot = "counts") +
  scale_fill_gradientn(colours = c('black', 'orange')) +
  scale_y_discrete(labels = genestegflip$protein)







### RECLUSTER TEGUMENT 


DefaultAssay(object = seurat_teg) <- "RNA"
seurat_teg <- DietSeurat(seurat_teg,counts = TRUE, data = FALSE ,scale.data = FALSE,features = NULL, assays = 'RNA', dimreducs = NULL, graphs = NULL)

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


DoHeatmap(object = seurat_teg, features = top20, group.by = "seurat_clusters", assay = 'RNA', slot = "scale.data") +
  scale_fill_gradientn(colours = c('black', 'orange')) +
  scale_y_discrete(labels = top20flipped$name)

#finds markers that disinguish the sub clusters from each other
allmarkersT <- FindAllMarkers(seurat_teg, assay = 'RNA', test.use =  'roc', only.pos = TRUE, return.thresh = 0)
markersT <- filter(allmarkersT, myAUC >= 0.8)

#find markers that distinguish the subclusters from every other cluster 




#write_xlsx(markersT,"markersT.xlsx")

markersannotatedT <- read.csv("markersTannotated.csv")
markersannotatedT <- markersannotated[,-10]
colnames(markersannotatedT) <- c('myAUC',"avg_diff","power","avg_logFC","pct.1","pct.2","cluster","gene","Gene.description")

markersannotatedT0 <- filter(markersannotatedT, cluster == 0)
markersannotatedT1 <- filter(markersannotatedT, cluster == 1)
markersannotatedT2 <- filter(markersannotatedT, cluster == 2)
markersannotatedT3 <- filter(markersannotatedT, cluster == 3)
markersannotatedT4 <- filter(markersannotatedT, cluster == 4)
markersannotatedT5 <- filter(markersannotatedT, cluster == 5)
markersannotatedT6 <- filter(markersannotatedT, cluster == 6)
markersannotatedT7 <- filter(markersannotatedT, cluster == 7)

DoHeatmap(object = seurat_teg, features = markersT$gene, group.by = "seurat_clusters", assay = 'RNA', slot = "scale.data") +
  scale_fill_gradientn(colours = c('blue', 'black', 'orange')) 

VlnPlot(seurat_object, features = 'Smp-340900', pt.size = 0.5 ) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
    labs(title = paste("NA","\n", "Smp-340900"))

#heatmap for tegument proteins on reclustered cells
DoHeatmap(object = seurat_teg, features = genesteg$gene_ID, group.by = "seurat_clusters", assay = 'RNA', slot = "scale.data") +
  scale_fill_gradientn(colours = c('blue', 'black', 'orange')) +
  scale_y_discrete(labels = genestegflip$protein)




