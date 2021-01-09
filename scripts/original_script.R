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
genes <- read.csv(file="gene list MASTER v2.csv")
names(genes) <- c("type", "protein", "gene_ID", "table", "tableName", "tableID", "stableID", "GeneName", "GeneDescription", "ignore")

#split the genes dataframe by which table in Alan's paper they came from
features1 <- subset(genes, genes$table==1)
features2 <- subset(genes, genes$table==2)
features3 <- subset(genes, genes$table==3)
#combine 2 and 3
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

library(ggplot2)
library(ggpubr)

####### plots for table 1
num <- 1: 10
for (i in num){
   #a <- features1$gene_ID[i]
   #png(sprintf("%s.png", a), width = 550, height = 350)
   assign(paste("Tone", i, sep = ''), FeaturePlot(seurat_object, features = features1$gene_ID[i]) +  
                labs(title = paste(features1$protein[i],"\n", features1$gene_ID[i]))+
                theme(plot.title = element_text(size=8)))
   #dev.off()
}

table1 <- ggarrange(Tone1, Tone2, Tone3, Tone4, Tone5, Tone6, Tone7, Tone8, Tone9, Tone10, ncol=3, nrow=3)
png("table1aUMAPs.png", width = 700, height = 700)
plot(table1$`1`)
dev.off()
png("table1bUMAPs.png", width = 700, height = 700)
plot(table1$`2`)
dev.off()


########## plots for table 2/3
num <- 1:66
for (i in num){
  #a <- features1$gene_ID[i]
  #png(sprintf("%s.png", a), width = 550, height = 350)
  assign(paste("Ttwo", i, sep = ''), FeaturePlot(seurat_object, features = features23$gene_ID[i]) +  
                                            labs(title = paste(features23$protein[i],"\n", features23$gene_ID[i]))+
                                            theme(plot.title = element_text(size=8)))
  #dev.off()
}

print(noquote(paste("Ttwo", 1:66, sep = '')))
table2 <- ggarrange(Ttwo1,  Ttwo2,  Ttwo3,  Ttwo4,  Ttwo5,  Ttwo6,  Ttwo7,  Ttwo8,  Ttwo9,  Ttwo10, 
                    Ttwo11, Ttwo12, Ttwo13, Ttwo14, Ttwo15, Ttwo16, Ttwo17, Ttwo18, Ttwo19, Ttwo20,
                    Ttwo21, Ttwo22, Ttwo23, Ttwo24, Ttwo25, Ttwo26, Ttwo27, Ttwo28, Ttwo29, Ttwo30,
                    Ttwo31, Ttwo32, Ttwo33, Ttwo34, Ttwo35, Ttwo36, Ttwo37, Ttwo38, Ttwo39, Ttwo40,
                    Ttwo41, Ttwo42, Ttwo43, Ttwo44, Ttwo45, Ttwo46, Ttwo47, Ttwo48, Ttwo49, Ttwo50,
                    Ttwo51, Ttwo52, Ttwo53, Ttwo54, Ttwo55, Ttwo56, Ttwo57, Ttwo58, Ttwo59, Ttwo60,
                    Ttwo61, Ttwo62, Ttwo63, Ttwo64, Ttwo65, Ttwo66, ncol=3, nrow=3)

png("table2aUMAPs.png", width = 700, height = 700)
plot(table2$`1`)
dev.off()
png("table2bUMAPs.png", width = 700, height = 700)
plot(table2$`2`)
dev.off()
png("table2cUMAPs.png", width = 700, height = 700)
plot(table2$`3`)
dev.off()
png("table2dUMAPs.png", width = 700, height = 700)
plot(table2$`4`)
dev.off()
png("table2eUMAPs.png", width = 700, height = 700)
plot(table2$`5`)
dev.off()
png("table2fUMAPs.png", width = 700, height = 700)
plot(table2$`6`)
dev.off()
png("table2gUMAPs.png", width = 700, height = 700)
plot(table2$`7`)
dev.off()
png("table2hUMAPs.png", width = 700, height = 700)
plot(table2$`8`)
dev.off()

######### plots for table 4
num <- 1: 24
for (i in num){
  #a <- features1$gene_ID[i]
  #png(sprintf("%s.png", a), width = 550, height = 350)
  assign(paste("Tfour", i, sep = ''), FeaturePlot(seurat_object, features = features4$gene_ID[i]) +  
                                            labs(title = paste(features4$protein[i],"\n", features4$gene_ID[i]))+
                                            theme(plot.title = element_text(size=8)))
  #dev.off()
}
print(noquote(paste("Tfour", 1:24, ",", sep = '')))
table4 <- ggarrange(Tfour1,  Tfour2,  Tfour3,  Tfour4,  Tfour5,  Tfour6,  Tfour7,  Tfour8, 
                    Tfour9,  Tfour10, Tfour11, Tfour12, Tfour13, Tfour14, Tfour15, Tfour16,
                    Tfour17, Tfour18, Tfour19, Tfour20, Tfour21, Tfour22, Tfour23, Tfour24,
                    ncol=3, nrow=3)

png("table4aUMAPs.png", width = 700, height = 700)
plot(table4$`1`)
dev.off()
png("table4bUMAPs.png", width = 700, height = 700)
plot(table4$`2`)
dev.off()
png("table4cUMAPs.png", width = 700, height = 700)
plot(table4$`3`)
dev.off()

######### plots for table 5
num <- 1: 8
for (i in num){
  #a <- features1$gene_ID[i]
  #png(sprintf("%s.png", a), width = 550, height = 350)
  assign(paste("Tfive", i, sep = ''), FeaturePlot(seurat_object, features = features5$gene_ID[i]) +  
           labs(title = paste(features5$protein[i],"\n", features5$gene_ID[i]))+
           theme(plot.title = element_text(size=8)))
  #dev.off()
}
print(noquote(paste("Tfive", 1:8, ",", sep = '')))
table5 <- ggarrange(Tfive1, Tfive2, Tfive3, Tfive4, Tfive5, 
                    Tfive6, Tfive7, Tfive8, nrow = 3, ncol =3)

png("table5aUMAPs.png", width = 700, height = 700)
plot(table5)
dev.off()

##END OF UMAPS

##### attempt to investigate transformation

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


### oesophageal gland work 

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
