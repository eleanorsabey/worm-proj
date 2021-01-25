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
top10
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
#feature plot for top 10 
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
data.frame(markers)

#put results of markers into wormbase and saved as top10 foldchange
#plot for genes with largest fold changes
FeaturePlot(object = seurat_object, cells = c, features = markers$rn, cols = c("orange", "blue"))

#Expression of three (MEGs 12, 16, 17) was confirmed in the anterior gland and 
#five (MEGs 8.1, 9, 11, 15 and 22) in the posterior gland.
#plot for genes shown to distinguish anterioir and posterior glands from literatire
eso_proteins <- read.csv("esophageal_gland_proteins.csv")
colnames(eso_proteins) <- c('name', 'ID')
num <- 1:71
for (i in num) {
  assign(paste('a', i, sep=''),FeaturePlot(object = seurat_object, cells = c, features = eso_proteins$ID[i], cols = c("orange", "blue")) +
           labs(title = paste(eso_proteins$name[i],"\n", eso_proteins$ID[i]))+
           theme(plot.title = element_text(size=8)))
}
print(noquote(paste("a", 1:74, ",", sep = '')))
eso <- ggarrange(a1,  a2,  a3,  a4,  a5,  a6,  a7,  a8,  a9,  a10, a11, a12, a13, a14,
                 a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28,
                 a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42,
                 a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56,
                 a57, a58, a59, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a70,
                 a71, a72, a73, a74, ncol=3, nrow=3)
png("esophageal_UMAPs.png", width = 700, height =700)
plot(eso$`1`)
dev.off()
png("esophageal_UMAPsb.png", width = 700, height =700)
plot(eso$`2`)
dev.off()
png("esophageal_UMAPsc.png", width = 700, height =700)
plot(eso$`3`)
dev.off()
png("esophageal_UMAPsd.png", width = 700, height =700)
plot(eso$`4`)
dev.off()
png("esophageal_UMAPse.png", width = 700, height =700)
plot(eso$`5`)
dev.off()
png("esophageal_UMAPsf.png", width = 700, height =700)
plot(eso$`6`)
dev.off()
png("esophageal_UMAPsg.png", width = 700, height =700)
plot(eso$`7`)
dev.off()
png("esophageal_UMAPsh.png", width = 700, height =700)
plot(eso$`8`)
dev.off()
png("esophageal_UMAPsi.png", width = 700, height =700)
plot(eso$`9`)
dev.off()

#extract coordinates of oesophageal cells
DimPlot(object = seurat_sub, cells = c, reduction = 'umap')

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

#add true neoblasts into a copy of the dataframe 
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
library(data.table)
setDT(coord2, keep.rownames = TRUE)[]
cluster <- c('neoblast')
neo_cells <- cbind(coord2, cluster)
colnames(neo_cells) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster')
neo_og_cells <- rbind(cluster_df, neo_cells)

### are the oesophageal cells not neoblasts
library(tidyr)
#change in Smp-ID for OG cell markers #meg8.2 = 'Smp-172180' this ID was used to highlight OG cells in the Wendt Paper
#megs 4.2 = Smp-085840 14 = Smp-124000 meant to mark OG gland  (https://doi.org/10.1371/journal.pntd.0002337)
#neoblast marker nanos2 = smp-051920
dat2 <- seurat_sub@assays[["integrated"]]@scale.data
dat3 <- seurat_sub_neo@assays[["integrated"]]@scale.data
dat2 <- data.frame(dat2)
dat3 <- data.frame(dat3)

newdat2 <- dat2['Smp-051920',]
newdat3 <- dat3['Smp-051920',]

cn2 <- colnames(newdat2)
cn3 <- colnames(newdat3)

setDT(newdat2, keep.rownames = TRUE)[]
setDT(newdat3, keep.rownames = TRUE)[]

ogcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
nbcountdf <- pivot_longer(newdat3, cols = cn3, names_to = "cell")

colnames(ogcountdf) <- c('ID', 'cell', 'scaled_counts_per_cell')
colnames(nbcountdf) <- c('ID', 'cell', 'scaled_counts_per_cell')

ognbcountdf <- rbind(ogcountdf, nbcountdf)

nc <- ognbcountdf$cell

ogandneo <- merge(ognbcountdf, neo_og_cells, by = 'cell')

nc <- neo_og_cells$cell

ggplot(data = ogandneo, aes(x = cluster, y = scaled_counts_per_cell)) +
  geom_violin() +
  geom_jitter(width=0.15, alpha=0.5)

#means for each cluster group
ogandneo %>% 
  group_by(cluster) %>% 
  summarise(mean(scaled_counts_per_cell))


nbcountsummary <-  ogandneo %>%
  group_by(cluster) %>%
  summarise(mean = mean(scaled_counts_per_cell),
            std = sd(scaled_counts_per_cell),
            n = length(scaled_counts_per_cell),
            se = std/sqrt(n))

#needs to be anova to compare the three clusters 
mod <- aov(scaled_counts_per_cell ~ cluster, data = ogandneo)
summary(mod)
TukeyHSD(mod)
plot(mod, 1)
plot(mod, 2)
#If data doesnt look normally distributed so use non-parametric test 

kw <- kruskal.test(scaled_counts_per_cell ~ cluster, data = ogandneo)
pairwise.wilcox.test(ogandneo$scaled_counts_per_cell, ogandneo$cluster,
                    p.adjust.method = "BH")



## heatmap for all og proteins 
c <- cluster_df$cell 
meta <- cluster_df[,-c(2:3)]
cell <- "ACTATCTCACGTTGGC_9"
cluster <- "outlier"
outlier <- data.frame(cell, cluster)
meta<- rbind(meta, outlier)
names <- meta[,1]
library(tidyr)
meta$cell <- as.factor(meta$cell)
meta$cluster <- as.factor(meta$cluster)
meta <- data.frame(meta)
row.names(meta) <- meta$cell
cluster.info <- meta 
seurat_sub <- AddMetaData(object = seurat_sub, metadata = cluster.info)

library(RColorBrewer)
#by default heatmap looks as scale.data slot, add cells = c to remove outlier column
eso_short <- eso_proteins
rownames(eso_short) <- eso_short$ID
eso_NA <- c('Smp-034420', 'Smp-132480', 'Smp-031180', 'Smp-300080', 'Smp-300070', 'Smp-167770', 'Smp-136240', 'Smp-213500', 'Smp-082030', 'Smp-071610', 'Smp-043390', 'Smp-019030', 'Smp-010620', 'Smp-002600', 'Smp-243780', 'Smp-243760', 'Smp-243750', 'Smp-331590', 'Smp-318200', 'Smp-266940', 'Smp-331380', 'Smp-165050', 'Smp-307240', 'Smp-307220', 'Smp-325680', 'Smp-138060', 'Smp-138070', 'Smp-138080', 'Smp-180340', 'Smp-180310', 'Smp-180320', 'Smp-159830', 'Smp-180330', 'Smp-159800')
esocleaned <- eso_short[!eso_short$ID %in% eso_NA, ]

#scale_y_discrete labels from opposite way to plotting so label from an inverted dataframe
library(purrr)
esoflipped <- esocleaned
esoflipped <- esoflipped %>% map_df(rev)


DoHeatmap(object = seurat_sub, features = esocleaned$ID, group.by = "cluster", draw.lines = TRUE)+ 
  scale_fill_gradientn(colours = c('blue', 'black', 'orange')) +
  scale_y_discrete(labels = esoflipped$name)
 

#try making heatmap read raw count expression slot

eso_proteinsflip <- eso_proteins
eso_proteinsflip <- eso_proteinsflip %>% map_df(rev)

DoHeatmap(object = seurat_sub, features = eso_proteins$ID, group.by = "cluster", assay = 'RNA', slot = "counts") +
  scale_fill_gradientn(colours = c('black', 'orange')) +
  scale_y_discrete(labels = eso_proteinsflip$name)

# 1. T test between expression of selected proteins between clusters of OG cells

#meg8.2 = 'Smp-172180' this ID was used to highlight OG cells in the Wendt Paper
#megs 4.2 = Smp-085840 14 = Smp-124000 meant to mark OG gland  (https://doi.org/10.1371/journal.pntd.0002337)

#Use scaled data slot
dat1 <- seurat_sub@assays[["integrated"]]@scale.data
dat1 <- data.frame(dat1)
#change Smp-ID for each 
newdat1 <- dat1['Smp-124000',]

cn1 <- colnames(newdat1)
setDT(newdat1, keep.rownames = TRUE)[]
scaledatadf <- pivot_longer(newdat1, cols = cn1, names_to = "cell")
colnames(scaledatadf) <- c('ID', 'cell', 'scaled_counts_per_cell')
mergedf1 <- merge(scaledatadf, cluster_df, by = "cell")

#t test to assess if theres a difference in ant meg expression
ggplot(data = mergedf1, aes(x = cluster, y = scaled_counts_per_cell)) +
  geom_violin()+
  geom_jitter(width=0.15, alpha=0.5)
#means for each cluster group
mergedf1 %>% 
  group_by(cluster) %>% 
  summarise(mean(scaled_counts_per_cell))


countsummary <- mergedf1 %>%
  group_by(cluster) %>%
  summarise(mean = mean(scaled_counts_per_cell),
            std = sd(scaled_counts_per_cell),
            n = length(scaled_counts_per_cell),
            se = std/sqrt(n))
t.test(data = mergedf1,
       scaled_counts_per_cell ~ cluster,
       var.equal = T)


# 2. Recluster the lower cluster and assess differences in expression of anterior/posterior MEGs.
# MEGs 12, 16, 17 should mark anterior gland 
#12 = Smp-152630, 16 = Smp-158890, 17 = Smp-180620


#recluster attempt

#lower cell df made earlier
colnames(lower) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster') 
lnames <- lower$cell

DimPlot(object = seurat_object, label = TRUE, pt.size = 0.5) + theme(legend.position = "none")

seurat_object[["ClusterNames_0.6"]] <- Idents(object = seurat_object)
seurat_object <- FindClusters(object = seurat_object, resolution = 1)
plot1 <- DimPlot(object = seurat_object, label = TRUE) + theme(legend.position = "none")
plot2 <- DimPlot(object = seurat_object, group.by = "ClusterNames_0.6", label = TRUE) + theme(legend.position = "none")
plot_grid(plot1, plot2)
plot3 <- DimPlot(object = seurat_object, label = TRUE, cells = nc)+geom_point(alpha=0.5)
plot4 <- DimPlot(object = seurat_object, group.by = "ClusterNames_0.6", label = TRUE, cells = nc)
plot4 + geom_point(alpha=0.5)
plot_grid(plot3, plot4)
DimPlot(object = seurat_object, reduction = 'umap', cells.highlight = c, label = TRUE) + theme(legend.position = "none") 
og.markers <- FindMarkers(object = seurat_object, ident.1 = 0, ident.2 = 54)
setDT(og.markers, keep.rownames = TRUE)[]

#doesn't work to recluster eso gland isolated 
#seurat_sub[["ClusterNames_0.6"]] <- Idents(object = seurat_sub)
#seurat_sub <- FindClusters(object = seurat_sub, resolution = 1)
#reluster attempt on og gland only ## desnt work
# DimPlot(object = seurat_sub, label = TRUE, pt.size = 0.5) + theme(legend.position = "none")
# seurat_sub[["ClusterNames_0.6"]] <- Idents(object = seurat_sub)
# seurat_sub <- FindClusters(object = seurat_sub, resolution = 0.8)
# plot1 <- DimPlot(object = seurat_sub, label = TRUE) + theme(legend.position = "none")
# plot2 <- DimPlot(object = seurat_sub, group.by = "ClusterNames_0.6", label = TRUE) + theme(legend.position = "none")
# plot_grid(plot1, plot2)

