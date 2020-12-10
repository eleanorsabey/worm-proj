library(Seurat)
#load data
seurat_object <- readRDS("data/GSE146736_adult_scseq_seurat.rds")
genes <- read.csv("gene list Master.csv")
colnames(genes) <- c('type','protein','gene_ID','table','tableName','tableID','genestableID','geneName','description','ignore')

DimPlot(seurat_object, reduction="umap")

dim(genes)
#set num to be the # of observations 
num <- 1:100

#make heatmaps, fontsize might not work, try + theme(plot.title = element_text(size=10))
for (i in num) {
  assign(paste('p', i, sep =''), FeaturePlot(seurat_object, features = genes$gene_ID[i]) +
  labs(title= paste(genes$protein[i],"\n", genes$gene_ID, sep = ''), size = 10))
}

library(ggplot2)
#remove key and rotate axis labels 90 degrees
for (i in num) {
  png(paste('v', i,'.png', sep = ''), width = 800, height = 500)
  plot(VlnPlot(seurat_object, features = genes$gene_ID[i]) +
    labs(title = paste(genes$protein[i], '/n', genes$gene_ID), size = 10))
  dev.off()
}
