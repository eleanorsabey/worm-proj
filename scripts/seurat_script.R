setwd("E:/project/scripts")
library(Seurat)
library(ggplot2)
library(ggpubr)

#load data
seurat_object <- readRDS("data/GSE146736_adult_scseq_seurat.rds")
all.genes <- rownames(seurat_object)
seurat_object<- ScaleData(seurat_object, features = all.genes)

DimPlot(seurat_object, reduction="umap", label = FALSE) + theme(legend.position = "none") 
genes <- read.csv("gene list Master v2.csv")
colnames(genes) <- c('type','protein','gene_ID','table','tableName','tableID','genestableID','geneName','description','ignore')

#whole seurat
DimPlot(seurat_object, reduction="umap")

max(seurat_object[["integrated"]]@scale.data)
min(seurat_object[["integrated"]]@scale.data)

hist(seurat_object@assays$integrated@scale.data)

#Making plots for indiviudal identifiers (vaccines in trial)
#Sm14 (Smp-095360,Smp-041430, should be in tegument and gut), 
#Sm-TSP-2 (Smp-335630 should be on tegument and likely in OG as this is extension of teg), 
#Calpain/Sm-p80 (Smp-214190 should be in tegument).

#530 by 430
FeaturePlot(seurat_object, features = 'Smp-041430',  slot = 'scale.data') + 
  scale_colour_gradientn(colours = c("#74acd1", "#dff3f8", "#fddf90", "#fcad60", "#F46D43", "#d62f27", "#a40025"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), na.value = 'white', name = paste('Expression', '\n', 'Z-score')) +
  labs(title = paste("Sm14 ","\n", "Smp-041430") )+
  theme(plot.title = element_text(size=10))

#980 by 400
VlnPlot(seurat_object, features = 'Smp-003300', pt.size = 0.5, slot = 'scale.data') +
       theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
       labs(title = paste("Sm14","\n", "Smp-041430"))
     
dim(genes)
#set num to be the # of observations 


#violin plots to produce indivual plot for each gene
num <- 1:108
for (i in num){
  png(sprintf("%s.png", paste('ScaleVln', i)), width = 850, height = 500)
  plot(VlnPlot(seurat_object, slot = 'scale.data', features = genes$gene_ID[i], pt.size = 0.5 ) +
         theme(legend.position = "none", axis.text.x = element_text(angle = 90))+
         labs(title = paste(genes$protein[i],"\n", genes$gene_ID[i])))
  dev.off()
}

num <- 1:108
for (i in num){
  png(sprintf("%s.png", paste('S', i)), width = 530, height = 430)
  plot(FeaturePlot(seurat_object, features = genes$gene_ID[i], slot = 'scale.data') + 
         scale_colour_gradientn(colours = c("#74acd1", "#dff3f8", "#fddf90", "#fcad60", "#F46D43", "#d62f27", "#a40025"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), na.value = 'white', name = paste('Expression', '\n', 'Z-score')) +
         labs(title = paste(genes$protein[i],"\n", genes$gene_ID[i]))+
    theme(plot.title = element_text(size=11)))
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
           scale_colour_gradientn(colours = c("#74acd1", "#dff3f8", "#fddf90", "#fcad60", "#F46D43", "#d62f27", "#a40025"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), na.value = 'white', name = paste('Expression', '\n', 'Z-score')) +
           labs(title = paste(features1$protein[i],"\n", features1$gene_ID[i]))+
           theme(plot.title = element_text(size=8)))
           
  
  #dev.off()
}

table1 <- ggarrange(Tone1, Tone2, Tone3, Tone4, Tone5, Tone6, Tone7, Tone8, Tone9, Tone10, ncol=3, nrow=3)
png("Stable1aUMAPs.png", width = 800, height = 700)
plot(table1$`1`)
dev.off()
png("Stable1bUMAPs.png", width = 800, height = 700)
plot(table1$`2`)
dev.off()


########## plots for table 2/3
num <- 1:66
for (i in num){
  #a <- features1$gene_ID[i]
  #png(sprintf("%s.png", a), width = 550, height = 350)
  assign(paste("Ttwo", i, sep = ''), FeaturePlot(seurat_object, features = features23$gene_ID[i]) +  
           scale_colour_gradientn(colours = c("#74acd1", "#dff3f8", "#fddf90", "#fcad60", "#F46D43", "#d62f27", "#a40025"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), na.value = 'white', name = paste('Expression', '\n', 'Z-score')) +
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

png("Stable2aUMAPs.png", width = 800, height = 700)
plot(table2$`1`)
dev.off()
png("Stable2bUMAPs.png", width = 800, height = 700)
plot(table2$`2`)
dev.off()
png("Stable2cUMAPs.png", width = 800, height = 700)
plot(table2$`3`)
dev.off()
png("Stable2dUMAPs.png", width = 800, height = 700)
plot(table2$`4`)
dev.off()
png("Stable2eUMAPs.png", width = 800, height = 700)
plot(table2$`5`)
dev.off()
png("Stable2fUMAPs.png", width = 800, height = 700)
plot(table2$`6`)
dev.off()
png("Stable2gUMAPs.png", width = 800, height = 700)
plot(table2$`7`)
dev.off()
png("Stable2hUMAPs.png", width = 800, height = 700)
plot(table2$`8`)
dev.off()

######### plots for table 4
num <- 1: 24
for (i in num){
  #a <- features1$gene_ID[i]
  #png(sprintf("%s.png", a), width = 550, height = 350)
  assign(paste("Tfour", i, sep = ''), FeaturePlot(seurat_object, features = features4$gene_ID[i]) +  
           scale_colour_gradientn(colours = c("#74acd1", "#dff3f8", "#fddf90", "#fcad60", "#F46D43", "#d62f27", "#a40025"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), na.value = 'white', name = paste('Expression', '\n', 'Z-score')) +
           labs(title = paste(features4$protein[i],"\n", features4$gene_ID[i]))+
           theme(plot.title = element_text(size=8)))
  #dev.off()
}
print(noquote(paste("Tfour", 1:24, ",", sep = '')))
table4 <- ggarrange(Tfour1,  Tfour2,  Tfour3,  Tfour4,  Tfour5,  Tfour6,  Tfour7,  Tfour8, 
                    Tfour9,  Tfour10, Tfour11, Tfour12, Tfour13, Tfour14, Tfour15, Tfour16,
                    Tfour17, Tfour18, Tfour19, Tfour20, Tfour21, Tfour22, Tfour23, Tfour24,
                    ncol=3, nrow=3)

png("Stable4aUMAPs.png", width = 800, height = 700)
plot(table4$`1`)
dev.off()
png("Stable4bUMAPs.png", width = 800, height = 700)
plot(table4$`2`)
dev.off()
png("Stable4cUMAPs.png", width = 800, height = 700)
plot(table4$`3`)
dev.off()

######### plots for table 5
num <- 1: 8
for (i in num){
  #a <- features1$gene_ID[i]
  #png(sprintf("%s.png", a), width = 550, height = 350)
  assign(paste("Tfive", i, sep = ''), FeaturePlot(seurat_object, features = features5$gene_ID[i]) +  
           scale_colour_gradientn(colours = c("#74acd1", "#dff3f8", "#fddf90", "#fcad60", "#F46D43", "#d62f27", "#a40025"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), na.value = 'white', name = paste('Expression', '\n', 'Z-score')) +
           labs(title = paste(features5$protein[i],"\n", features5$gene_ID[i]))+
           theme(plot.title = element_text(size=8)))
  #dev.off()
}
print(noquote(paste("Tfive", 1:8, ",", sep = '')))
table5 <- ggarrange(Tfive1, Tfive2, Tfive3, Tfive4, Tfive5, 
                    Tfive6, Tfive7, Tfive8, nrow = 3, ncol =3)

png("Stable5aUMAPs.png", width = 800, height = 700)
plot(table5)
dev.off()

##END OF UMAPS


teg <- seurat_object@active.ident
tegdf <- data.frame(teg)
setDT(tegdf, keep.rownames = TRUE)[]
colnames(tegdf) <- c("cell", "tissue")
tegdf$tissue <- as.factor(tegdf$tissue)
teg1 <- filter(tegdf, tegdf$tissue == "tegument 1")
teg2 <- filter(tegdf, tegdf$tissue == "tegument 2")
teg3 <- filter(tegdf, tegdf$tissue == "early tsp-2+")
teg4 <- filter(tegdf, tegdf$tissue == "egc+")
teg5 <- filter(tegdf, tegdf$tissue == "meg-1+")
teg6 <- filter(tegdf, tegdf$tissue == "zfp-1-1+")
teg7 <- filter(tegdf, tegdf$tissue == "sm13+")

teg <- rbind(teg1, teg2, teg3, teg4, teg5, teg6, teg7)
c <- teg$cell



#subset the larger seurat object into a smaller one just for tegument
seurat_allteg <- subset(seurat_object, cells = c)

# all tegument progenitor cells ("early tsp-2+", "egc+", "meg-1+", "zfp-1-1+", and "sm13+") into a



#separate into upper and lower clusters, label, then combine
library(dplyr)


cluster_df <- rbind(teg1, teg2, teg3, teg4, teg5, teg6, teg7)
colnames(cluster_df) <- c('cell', 'tissue')
cluster <- c('tegument 1', 'tegument 2', 'early tsp-2+', "egc+", "meg-1+", "zfp-1-1+", "sm13+")

dat2 <- seurat_allteg@assays[["integrated"]]@scale.data
dat2 <- data.frame(dat2)

#genes to be interested in:
#Making plots for indiviudal identifiers (vaccines in trial)
#Sm14 (Smp-095360,Smp-041430, should be in tegument and gut), 
#Sm-TSP-2 (Smp-335630 should be on tegument and likely in OG as this is extension of teg), 
#Calpain/Sm-p80 (Smp-214190 should be in tegument).
genes <- c('Smp-041430','Smp-335630','Smp-214190')
names <- c('Sm14','Tsp2','Calpain')
stat <- cbind(genes, names)
stat <- data.frame(stat)
num <- 1:3



names <- list()
ids <- list()
teg1 <- list()
teg2 <- list()
earlytsp2 <- list()
egc <- list()
meg1 <- list()
sm13 <- list()
zfp_1_1 <- list()

#use stat for vaccine canididates, features23 for teg protein, features4 for vomitus
for (i in num){
  
  n <- stat$names[i]
  names <- append(names, n)
  id <- stat$genes[i]
  ids <- append(ids, id)
  
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  tegcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(tegcountdf) <- c('ID', 'cell', 'expression')
  
  tegandprg <- merge(tegcountdf, cluster_df, by = 'cell')
  
  #Mean expression for each group:
  #print('Mean expression:')
  m <- tegandprg %>% 
    group_by(tissue) %>% 
    summarise(mean(expression))
  
  t2 <- m$`mean(expression)`[1]
  teg2 <- append(teg2, t2)
  t1 <- m$`mean(expression)`[2]
  teg1 <- append(teg1, t1)
  sm <- m$`mean(expression)`[3]
  sm13 <- append(sm13,sm)
  e <- m$`mean(expression)`[4]
  egc <- append(egc,e)
  t <- m$`mean(expression)`[5]
  earlytsp2 <- append(earlytsp2,t)
  meg <- m$`mean(expression)`[6]
  meg1 <- append(meg1,meg)
  z <- m$`mean(expression)`[7]
  zfp_1_1 <- append(zfp_1_1,z)
}
num <- 1:66
for (i in num){
  
  n <- features23$protein[i]
  names <- append(names, n)
  id <- features23$gene_ID[i]
  ids <- append(ids, id)
  
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  tegcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(tegcountdf) <- c('ID', 'cell', 'expression')
  
  tegandprg <- merge(tegcountdf, cluster_df, by = 'cell')
  
  #Mean expression for each group:
  #print('Mean expression:')
  m <- tegandprg %>% 
    group_by(tissue) %>% 
    summarise(mean(expression))
  
  t2 <- m$`mean(expression)`[1]
  teg2 <- append(teg2, t2)
  t1 <- m$`mean(expression)`[2]
  teg1 <- append(teg1, t1)
  sm <- m$`mean(expression)`[3]
  sm13 <- append(sm13,sm)
  e <- m$`mean(expression)`[4]
  egc <- append(egc,e)
  t <- m$`mean(expression)`[5]
  earlytsp2 <- append(earlytsp2,t)
  meg <- m$`mean(expression)`[6]
  meg1 <- append(meg1,meg)
  z <- m$`mean(expression)`[7]
  zfp_1_1 <- append(zfp_1_1,z)
}

s_df <- rbind(names,
              ids,
              teg1,
              teg2,
              sm13,
              egc,
              meg1,
              zfp_1_1
              )

s_df <- data.frame(s_df)
s_df<- rotate_df(s_df)

#OESOPHAGEAL GLAND
og <- filter(tegdf, tegdf$tissue == "oesophageal gland")
gut <- filter(tegdf, tegdf$tissue == "gut")
og_gut <- rbind(og, gut)
c <- og_gut$cell
seurat_og_gut <- subset(seurat_object, cells = c)
cluster_df <- rbind(og, gut)
colnames(cluster_df) <- c('cell', 'tissue')
cluster <- c('oesophageal gland', 'gut')
dat2 <- seurat_og_gut@assays[["integrated"]]@scale.data
dat2 <- data.frame(dat2)
names <- list()
ids <- list()
oesophageal_gland <- list()
gut <- list()

num <- 1:24
for (i in num){
  
  n <- features4$protein[i]
  names <- append(names, n)
  id <- features4$gene_ID[i]
  ids <- append(ids, id)
  
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  ogcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(ogcountdf) <- c('ID', 'cell', 'expression')
  
  ogandgut <- merge(ogcountdf, cluster_df, by = 'cell')
  
  #Mean expression for each group:
  #print('Mean expression:')
  m <- ogandgut %>% 
    group_by(tissue) %>% 
    summarise(mean(expression))
  
  og <- m$`mean(expression)`[1]
  oesophageal_gland <- append(oesophageal_gland, og)
  g <- m$`mean(expression)`[2]
  gut <- append(gut, g)
  
}

s_df <- rbind(names,
              ids,
              oesophageal_gland,
              gut
)

s_df <- data.frame(s_df)
s_df<- rotate_df(s_df)
