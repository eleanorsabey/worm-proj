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

#doesn't work to recluster eso gland isolated 
#seurat_sub[["ClusterNames_0.6"]] <- Idents(object = seurat_sub)
#seurat_sub <- FindClusters(object = seurat_sub, resolution = 1)


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
num <- 1:74
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

#remove single cell outlier
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

# MEGs 12, 16, 17 should mark anterior gland 
#12 = Smp-152630, 16 = Smp-158890, 17 = Smp-180620
library(tidyr)
antMEGS <- c('Smp-152630', 'Smp-158890', 'Smp-180620')
#raw count data
dat <- seurat_sub@assays[["RNA"]]@counts
dat <- data.frame(dat)
newdat <- dat[antMEGS,]
cn <- colnames(newdat)
setDT(newdat, keep.rownames = TRUE)[]
countdf <- pivot_longer(newdat, cols = cn, names_to = "cell")
colnames(countdf) <- c('ID', 'cell', 'rna_per_cell')
mergedf <- merge(countdf, cluster_df, by = "cell")

#scaled data
dat1 <- seurat_sub@assays[["integrated"]]@scale.data
dat1 <- data.frame(dat1)
newdat1 <- dat1[antMEGS,]
cn1 <- colnames(newdat1)
setDT(newdat1, keep.rownames = TRUE)[]
scaledatadf <- pivot_longer(newdat1, cols = cn1, names_to = "cell")
colnames(scaledatadf) <- c('ID', 'cell', 'scaled_counts_per_cell')
mergedf1 <- merge(scaledatadf, cluster_df, by = "cell")

#t test to assess if theres a difference in ant meg expression
ggplot(data = mergedf, aes(x = cluster, y = rna_per_cell)) +
  geom_violin()
#means for each cluster group
mergedf %>% 
  group_by(cluster) %>% 
  summarise(mean(rna_per_cell))
#1 lower                  0.232
#2 upper                  0.222

countsummary <- mergedf %>%
  group_by(cluster) %>%
  summarise(mean = mean(rna_per_cell),
            std = sd(rna_per_cell),
            n = length(rna_per_cell),
            se = std/sqrt(n))
t.test(data = mergedf,
       rna_per_cell ~ cluster,
       var.equal = T)
#no significant difference between expression of anteriir meg
#RAW: t = 0.060608, df = 172, p-value = 0.9517


ggplot(data = mergedf1, aes(x = cluster, y = scaled_counts_per_cell)) +
  geom_violin()
#means for each cluster group
mergedf1 %>% 
  group_by(cluster) %>% 
  summarise(mean(scaled_counts_per_cell))
#1 lower                  -0.0471
#2 upper                  -0.0454

countsummary <- mergedf1 %>%
  group_by(cluster) %>%
  summarise(mean = mean(scaled_counts_per_cell),
            std = sd(scaled_counts_per_cell),
            n = length(scaled_counts_per_cell),
            se = std/sqrt(n))
t.test(data = mergedf1,
       scaled_counts_per_cell ~ cluster,
       var.equal = T)
# no significant difference 
#SCALED: t = -0.015033, df = 172, p-value = 0.988


#are the oesophageal cells not neoblasts
#neoblast marker nanos2 = smp-051920
dat2 <- seurat_sub@assays[["RNA"]]@counts
dat2 <- data.frame(dat2)
newdat2 <- dat2['Smp-051920',]
cn2 <- colnames(newdat2)
setDT(newdat2, keep.rownames = TRUE)[]
nbcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
colnames(nbcountdf) <- c('ID', 'cell', 'counts_per_cell')
mergedf2 <- merge(nbcountdf, cluster_df, by = "cell")

ggplot(data = mergedf2, aes(x = cluster, y = counts_per_cell)) +
  geom_violin()
#means for each cluster group
mergedf2 %>% 
  group_by(cluster) %>% 
  summarise(mean(counts_per_cell))


nbcountsummary <- mergedf2 %>%
  group_by(cluster) %>%
  summarise(mean = mean(counts_per_cell),
            std = sd(counts_per_cell),
            n = length(counts_per_cell),
            se = std/sqrt(n))
t.test(data = mergedf2,
       counts_per_cell ~ cluster,
       var.equal = T)
#RAW: t = -2.3389, df = 56, p-value = 0.02293
#SCALED: t = -3.0284, df = 56, p-value = 0.003714

#difference in expression for all eso proteins between upper and lower 
library(tidyr)
eso <- eso_proteins$ID
#make dataframe 
dat <- seurat_sub@assays[["integrated"]]@scale.data
dat <- data.frame(dat)
newdat <- dat[eso,]
cn <- colnames(newdat)
setDT(newdat, keep.rownames = TRUE)[]
countdf <- pivot_longer(newdat, cols = cn, names_to = "cell")
colnames(countdf) <- c('ID', 'cell', 'rna_per_cell')
mergedf <- merge(countdf, cluster_df, by = "cell")
#analyse
ggplot(data = mergedf, aes(x = cluster, y = rna_per_cell)) +
  geom_violin()
#means for each cluster group
mergedf %>% 
  group_by(cluster) %>% 
  summarise(mean(rna_per_cell))
#1 lower                  R= 0.718, S= 2.0094990
#2 upper                  R= 0.187, S= 0.4549808 

countsummary <- mergedf %>%
  group_by(cluster) %>%
  summarise(mean = mean(rna_counts_per_cell),
            std = sd(rna_counts_per_cell),
            n = length(rna_counts_per_cell),
            se = std/sqrt(n))

t.test(data = mergedf,
       rna_per_cell ~ cluster,
       var.equal = T)
#RAW:t = 5.179, df = 4290, p-value = 2.333e-07
#SCALED:t = 7.9612, df = 2202, p-value = 2.703e-15

#reluster attempt on og gland only ## desnt work
# DimPlot(object = seurat_sub, label = TRUE, pt.size = 0.5) + theme(legend.position = "none")
# seurat_sub[["ClusterNames_0.6"]] <- Idents(object = seurat_sub)
# seurat_sub <- FindClusters(object = seurat_sub, resolution = 0.8)
# plot1 <- DimPlot(object = seurat_sub, label = TRUE) + theme(legend.position = "none")
# plot2 <- DimPlot(object = seurat_sub, group.by = "ClusterNames_0.6", label = TRUE) + theme(legend.position = "none")
# plot_grid(plot1, plot2)

