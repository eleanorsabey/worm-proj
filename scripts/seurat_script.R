library(Seurat)
library(ggplot2)
library(ggpubr)

#load data
seurat_object <- readRDS("data/GSE146736_adult_scseq_seurat.rds")

DimPlot(seurat_object, reduction="umap", label = FALSE) + theme(legend.position = "none") 
genes <- read.csv("gene list Master v2.csv")
colnames(genes) <- c('type','protein','gene_ID','table','tableName','tableID','genestableID','geneName','description','ignore')

#whole seurat
DimPlot(seurat_object, reduction="umap")

max(seurat_object[["integrated"]]@scale.data)
min(seurat_object[["integrated"]]@scale.data)


hist(seurat_object@assays$integrated@scale.data)


FeaturePlot(seurat_object, features = 'Smp-105360',  slot = 'scale.data') + 
  scale_colour_gradientn(colours = c("#ABD9E9", "#E0F3F8", "#FDAE61", "#F46D43", "#D73027"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), name = paste('Expression', '\n', 'Z-score')) +
  labs(title = paste("notch","\n", "Smp-105360") )+
  theme(plot.title = element_text(size=10))

VlnPlot(seurat_object, features = 'Smp-090080', pt.size = 0.5 ) +
       theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
       labs(title = paste("Serpin","\n", "Smp-090080"))
     
dim(genes)
#set num to be the # of observations 


#violin plots to produce indivual plot for each gene
num <- 1:108
for (i in num){
  png(sprintf("%s.png", paste('v', i)), width = 850, height = 500)
  plot(VlnPlot(seurat_object, features = genes$gene_ID[i], pt.size = 0.5 ) +
         theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
         labs(title = paste(genes$protein[i],"\n", genes$gene_ID[i])))
  dev.off()
}



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



####### plots for table 1
num <- 1: 10
for (i in num){
  #a <- features1$gene_ID[i]
  #png(sprintf("%s.png", a), width = 550, height = 350)
  assign(paste("Tone", i, sep = ''), FeaturePlot(seurat_object, features = features1$gene_ID[i]) +  
           scale_colour_gradientn(colours = c("#ABD9E9", "#E0F3F8", "#FDAE61", "#F46D43", "#D73027"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10)) +
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
           scale_colour_gradientn(colours = c("#ABD9E9", "#E0F3F8", "#FDAE61", "#F46D43", "#D73027"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10))) +
           labs(title = paste(features23$protein[i],"\n", features23$gene_ID[i]))+
           theme(plot.title = element_text(size=8)))
  #dev.off()
}

#print(noquote(paste("Ttwo", 1:66, sep = '')))
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
           scale_colour_gradientn(colours = c("#ABD9E9", "#E0F3F8", "#FDAE61", "#F46D43", "#D73027"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), name = paste('Expression', '\n', 'Z-score')) +
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
           scale_colour_gradientn(colours = c("#ABD9E9", "#E0F3F8", "#FDAE61", "#F46D43", "#D73027"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), name = paste('Expression', '\n', 'Z-score')) +
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