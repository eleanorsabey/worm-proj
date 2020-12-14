setwd("E:/project/scripts")
invisible(utils::memory.limit(50000))
# Load Seurat into R
library(Seurat)

# Load object
seurat_object <- readRDS("data/GSE146736_adult_scseq_seurat.rds")

# Have a look at the object
seurat_object

# Plot the UMAP projection to see something very similar to Figure 1A of the paper
DimPlot(seurat_object, reduction="umap")



#load gene of interest list table in and adjust col names
genes <- read.csv(file="gene list MASTER.csv")
names(genes) <- c("type", "protein", "gene_ID", "table", "tableName", "tableID", "stableID", "GeneName", "GeneDescription", "ignore")

#split the genes dataframe by which table in Alan's paper they came from
features1 <- subset(genes, genes$table==1)
features2 <- subset(genes, genes$table==2)
features3 <- subset(genes, genes$table==3)
features23 <- rbind(features2, features3)
features4 <- subset(genes, genes$table==4)
features5 <- subset(genes, genes$table==5)
#split 23 and 4 into smaller frames to make plots look nicer
features4a <- features4[1:9,]
features4b <- features4[10:18,]
features4c <- features4[19:23,]
features23a <- features23[1:9,]
features23b <- features23[10:18,]
features23c <- features23[19:27,]
features23d <- features23[28:36,]
features23e <- features23[37:45,]
features23f <- features23[46:54,]
features23g <- features23[55:63,]

#make feature plots
library(ggplot2)
png("table 1.png", width = 700, height = 700)
plot(FeaturePlot(seurat_object, features = features1$gene_ID) +
       labs(title = paste(features1$protein, "\n",features1$gene_ID)))
dev.off()
DefaultDimReduc(object =features1$gene_ID )
num <- 1:100
for (n in num){
  
}

n14 <- FeaturePlot(object = seurat_object, features = 'Smp-313560') +
  labs(title = 'Alkaline phosphatase Smp-313560') +
  theme(plot.title = element_text(size=15))
plot(n14)
#violin plots to produce indivual plot for each gene
for(g in genes$gene_ID){
      png(sprintf("%s.png", , width = 850, height = 500)  
      plot(VlnPlot(seurat_object, features = g, pt.size = 0.5 ) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90)))
      dev.off()
}

for (i in num){
  #a <- features1$gene_ID[i]
  #png(sprintf("%s.png", a), width = 550, height = 350)
  assign(paste("v", i, sep = ''), plot(VlnPlot(seurat_object, features = genes$gene_ID[i], pt.size = 0.5 ) +
                                         theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
                                         labs(title = paste(genes$protein[i],"\n", genes$gene_ID[i]))))
 
  #dev.off()
}

for (i in num){
  png(sprintf("%s.png", paste('v', i)), width = 850, height = 500)
  plot(VlnPlot(seurat_object, features = genes$gene_ID[i], pt.size = 0.5 ) +
         theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
         labs(title = paste(genes$protein[i],"\n", genes$gene_ID[i])))
  dev.off()
}

plot(v1)
#loop for heatmap plots
#  num <- 1: 100 
# for (i in num){
#      #a <- genes$gene_ID[i]
#      #png(sprintf("%s.png", a), width = 550, height = 350)
#      plot(FeaturePlot(seurat_object, features = genes$gene_ID[i]) +  
#      labs(title = paste(genes$protein[i],genes$gene_ID[i])))
#      #dev.off()
# }

 num <- 1: 100


for (i in num){
   #a <- features1$gene_ID[i]
   #png(sprintf("%s.png", a), width = 550, height = 350)
   assign(paste("p", i, sep = ''), plot(FeaturePlot(seurat_object, features = genes$gene_ID[i]) +  
                labs(title = paste(genes$protein[i],"\n", genes$gene_ID[i]))+
                theme(plot.title = element_text(size=8))))
   #dev.off()
}

 table2anew <- ggarrange(p8,  p9,  p10, p11, p12, p13, n14, p15, p16, ncol = 3, nrow = 3)
 
library(ggpubr)
table1 <- ggarrange(p1, p2, p3, p4, p5, p6, p7,ncol = 3, nrow = 3) 
table2a <- ggarrange(p8,  p9,  p10, p11, p12, p13, p14, p15, p16, ncol = 3, nrow = 3)
table2b <- ggarrange(p17, p18, p19, p20, p21, p22, p23, p24, p25, ncol = 3, nrow = 3)
table2c <- ggarrange(p26, p28, p29, p30, p31, p32, p33, p34, p35, ncol = 3, nrow =3)
table2d <- ggarrange(p36, p38, p39, p40, p41, p42, p43, p44, p45,ncol = 3, nrow =3)
table2e <- ggarrange(p46, p51, p52, p53, p54, p55, p56, p57, p58,ncol = 3, nrow =3)
table2f <- ggarrange(p59, p60, p61, p62, p27, p37, p47, p48, p49,ncol = 3, nrow =3)
table2g <- ggarrange(p50, p63, p64, p65, p66, p67, p68, p69,ncol = 3, nrow =3)
table4a <- ggarrange(p70, p71, p72, p73, p74, p75, p76, p77, p78,ncol = 3, nrow =3)
table4b <- ggarrange(p79, p80, p81, p82, p83, p84, p85, p86, p87,ncol = 3, nrow =3)
table4c <- ggarrange(p88, p89, p90, p91, p92,ncol = 3, nrow =3)
table5 <- ggarrange(p93,  p94,  p95,  p96,  p97,  p98,  p99,  p100,ncol = 3, nrow =3)

png("table 1.png", width = 700, height = 700)
plot(table1)
dev.off()
png("table 2a.png", width = 700, height = 700)
plot(table2a)
dev.off()
png("table 2b.png", width = 700, height = 700)
plot(table2b)
dev.off()
png("table 2c.png", width = 700, height = 700)
plot(table2c)
dev.off()
png("table 2d.png", width = 700, height = 700)
plot(table2d)
dev.off()
png("table 2e.png", width = 700, height = 700)
plot(table2e)
dev.off()
png("table 2f.png", width = 700, height = 700)
plot(table2f)
dev.off()
png("table 2g.png", width = 700, height = 700)
plot(table2g)
dev.off()
png("table 4a.png", width = 700, height = 700)
plot(table4a)
dev.off()
png("table 4b.png", width = 700, height = 700)
plot(table4b)
dev.off()
png("table 4c.png", width = 700, height = 700)
plot(table4c)
dev.off()
png("table 5.png", width = 700, height = 700)
plot(table5)
dev.off()




list(features1$gene_ID)
#distributions
hist(seurat_object@assays$integrated@scale.data)
# > max(seurat_object@assays$integrated@scale.data)
# [1] 10
# > min(seurat_object@assays$integrated@scale.data)
# [1] -26.72586
# > 2^10
# [1] 1024
# > 2^-26.72586
# [1] 9.009786e-09

hist(seurat_object@assays$integrated@scale.data['Smp-051920',])
hist(seurat_object@assays$integrated@data@p)
hist(seurat_object@assays$integrated['Smp-051920',])
hist(seurat_object@assays$integrated@data@i)
hist(seurat_object@reductions$umap@cell.embeddings)
hist(seurat_object@meta.data$nCount_RNA)
hist(seurat_object@assays$RNA@counts@p)
hist(seurat_object@assays$integrated@data@x)
FeaturePlot(seurat_object, features = 'Smp-051920')
max(seurat_object@assays$integrated@scale.data['Smp-051920',])
#3.503137
min(seurat_object@assays$integrated@scale.data['Smp-051920',])
#-2.602912
min(seurat_object@assays$integrated['Smp-051920',])
#-2.080313
min(seurat_object@assays$integrated@data@x)
#-3.944866
max(seurat_object@assays$integrated@data@x)
#8.841464
class(seurat_object@assays$integrated@data@x)

hist(seurat_object@assays$RNA@counts@x)
min(seurat_object@assays$RNA@counts@x)
#1
max(seurat_object@assays$RNA@counts@x)
#8007
mean(seurat_object@assays$RNA@counts@x)
#2.413571

orig <- seurat_object@assays$RNA@counts@x
hist(orig)

logged <- log2(orig)
hist(logged)

final <- seurat_object@assays$integrated@data@x
hist(final)
hist(seurat_object@assays$integrated@data@x)
min(seurat_object@assays$integrated@data@x)
#-3.944866
max(seurat_object@assays$integrated@data@x)
#8.841464
mean(seurat_object@assays$integrated@data@x)
#0.3868306

library(dplyr)
library(Seurat)
library(patchwork)

normalised <- seurat_object[["RNA"]]@data
raw0 <- seurat_object[["RNA"]]@counts
raw <- NormalizeData(raw0)

final <- seurat_object[["integrated"]]@data

raw1 <- NormalizeData(raw)
og <- seurat_object@active.ident
DimPlot(object = seurat_object, cells = c, reduction = 'umap')
FeaturePlot(object = seurat_object, cells = c, features = 'Smp-067060')
?FeaturePlot
?DimPlot
ogdf <- data.frame(og)
library(data.table)
setDT(ogdf, keep.rownames = TRUE)[]
colnames(ogdf) <- c("cell", "tissue")
ogdf$tissue <- as.factor(ogdf$tissue)
library(dplyr)
d <- filter(ogdf, ogdf$tissue == "oesophageal gland")
c <- d$cell

MEGs <- read.csv("MEGs.csv")
colnames(MEGs) <- c("name", "ID")
num0 <- 1:37
library(ggplot2)
for (i in num0){
  #a <- MEGs$ID[i]
  png(paste("m", i, ".png", sep = ''), width = 550, height = 350)
  plot(FeaturePlot(object = seurat_object, cells = c, features = MEGs$ID[i]) +  
                                         labs(title = paste(MEGs$name[i],"\n", MEGs$ID[i]))+
                                         theme(plot.title = element_text(size=10)))
  dev.off()
}

rm(d)
rm(og)
rm(ogdf)

#extract top 10 identifiers for OG cluster
cluster1.markers <- FindMarkers(seurat_object, ident.1 = c, min.pct = 0.25)
head(cluster1.markers, n = 10)
int <- seurat_object@assays[["integrated"]]
rm(seurat_object)

so1 <- SubsetData(object = seurat_object, assay = "integrated")
pbmc1
