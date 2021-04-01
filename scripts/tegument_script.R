#Tegument analysis

#####
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
library(sjmisc)

#####
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

#subset the larger seurat object into a smaller one just for tegument
seurat_teg <- subset(seurat_object, cells = c)
all.genes <- rownames(seurat_teg)
seurat_teg <- ScaleData(seurat_teg, features = all.genes)
# all tegument progenitor cells ("early tsp-2+", "egc+", "meg-1+", "zfp-1-1+", and "sm13+") into a
# single cluster ("Tegument Progenitors"),




#map of all ce;;s
DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + theme(legend.position = "none") 

#highlight cells labelled as tegument amongst all cells
DimPlot(object = seurat_object, reduction = 'umap', cells.highlight = c) + theme(legend.position = "none") 
#plot only tegument cells, this shows how the two tegument identities locate
DimPlot(object = seurat_object, cells = c, reduction = 'umap')
 

#####
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
seurat_neo <- subset(seurat_object, cells = c) 

coord <- seurat_neo[["umap"]]@cell.embeddings
coord2 <- as.data.frame(coord)
setDT(coord2, keep.rownames = TRUE)[]
cluster <- c('neoblast')
neo_cells <- cbind(coord2, cluster)
colnames(neo_cells) <- c('cell', 'UMAP_1', 'UMAP_2', 'cluster')
neo_teg_cells <- rbind(cluster_df, neo_cells)

#make a seurat object with teg and Neo cells, then re scale.
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
num <- 1:5

fileConn<-file("teg_vs_neo_stats.txt")
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
  assign(paste("hist", i, sep = ''), as.ggplot(~qqPlot(mod$residuals)))
  #Check equal varience:
  assign(paste("box", i, sep = ''), as.ggplot(~boxplot(expression ~ cluster, data = tegandneo)))
  
  
  print("Non-parametric - Kruskal Wallis and unpaired Wilcox test:")
  #If varience isnt equal: Welch's anova 
  print(kruskal.test(expression ~ cluster, data = tegandneo))
  print(pairwise.wilcox.test(tegandneo$expression, tegandneo$cluster,
                             p.adjust.method = "none", paired = FALSE))
  
}
sink()
close(fileConn)


statplot <- ggarrange(Plot1, hist1, box1, Plot2, hist2, box2, Plot3, hist3, box3, Plot4,  hist4, box4, ncol=3, nrow=4)
statplot2 <- ggarrange(Plot5, hist5, box5, ncol=3, nrow=4)
plot(statplot)
png("teg_vs_neo_stats_1.png", width = 470, height = 800)
plot(statplot)
dev.off()
png("teg_vs_neo_stats_2.png", width = 470, height = 800)
plot(statplot2)
dev.off()
plot(Plot1)

##### write stats to table
names <- list()
ids <- list()
lower1 <- list()
neoblast <- list()
upper1 <- list()
upper2 <- list()
FvalANOVA <- list()
PvalANOVA <- list()
Tukey_Pval_neoblast_lower1 <- list()
Tukey_Pval_upper1_lower1 <- list()
Tukey_Pval_upper2_lower1 <- list()
Tukey_Pval_upper1_neoblast <- list()
Tukey_Pval_upper2_neoblast <- list()
Tukey_Pval_upper2_upper1 <- list()
KruskalWallis <- list()
PvalKruskal <- list()
Wilcox_Pval_neoblast_lower1 <- list()
Wilcox_Pval_upper1_lower1 <- list()
Wilcox_Pval_upper2_lower1 <- list()
Wilcox_Pval_upper1_neoblast <- list()
Wilcox_Pval_upper2_neoblast <- list()
Wilcox_Pval_upper2_upper1 <- list()

for (i in num){
  
  n <- stat$names[i]
  names <- append(names, n)
  id <- stat$genes[i]
  ids <- append(ids, id)
  

  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  tegnbcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(tegnbcountdf) <- c('ID', 'cell', 'expression')
  
  tegandneo <- merge(tegnbcountdf, neo_teg_cells, by = 'cell')
  
  #Mean expression for each group:
  #print('Mean expression:')
  m <- tegandneo %>% 
          group_by(cluster) %>% 
          summarise(mean(expression))
  
  l1 <- m$`mean(expression)`[1]
  lower1 <- append(lower1, l1)
  n <- m$`mean(expression)`[2]
  neoblast <- append(neoblast, n)
  u1 <- m$`mean(expression)`[3]
  upper1 <- append(upper1,u1)
  u2 <- m$`mean(expression)`[4]
  upper2 <- append(upper2,u2)
  
  #Anova to compare expression across the three clusters:
  mod <- aov(expression ~ cluster, data = tegandneo)
  #View anova:
  sum <- summary(mod)
  f <- sum[[1]][["F value"]][1]
  FvalANOVA <- append(f, FvalANOVA)
  p <- sum[[1]][["Pr(>F)"]][1]
  PvalANOVA <- append(p, PvalANOVA)
  
  tukey <- TukeyHSD(mod)
  tukeydf <- data.frame(tukey$cluster)
  n_l1 <- tukeydf$p.adj[1]
  Tukey_Pval_neoblast_lower1<- append(Tukey_Pval_neoblast_lower1, n_l1)
  
  u1_l1 <- tukeydf$p.adj[2]
  Tukey_Pval_upper1_lower1 <- append(Tukey_Pval_upper1_lower1, u1_l1)
  
  u2_l1 <- tukeydf$p.adj[2]
  Tukey_Pval_upper2_lower1 <- append(Tukey_Pval_upper2_lower1, u2_l1)
  
  u1_n <- tukeydf$p.adj[3]
  Tukey_Pval_upper1_neoblast <- append(Tukey_Pval_upper1_neoblast, u1_n)
  
  u2_n <- tukeydf$p.adj[3]
  Tukey_Pval_upper2_neoblast <- append(Tukey_Pval_upper2_neoblast, u2_n)
  
  u2_u1 <- tukeydf$p.adj[2]
  Tukey_Pval_upper2_upper1 <- append(Tukey_Pval_upper2_upper1, u2_u1)
  
  k<- kruskal.test(expression ~ cluster, data = tegandneo)
  chi2 <- k[["statistic"]][["Kruskal-Wallis chi-squared"]]
  KruskalWallis <- append(KruskalWallis, chi2)
  pkw <- k[["p.value"]]
  PvalKruskal <- append(PvalKruskal, pkw)
  
  w <- pairwise.wilcox.test(tegandneo$expression, tegandneo$cluster,
                             p.adjust.method = "none", paired = FALSE)
  wdf <- data.frame(w$p.value)
  w_nl1 <- wdf$lower1[1]
  Wilcox_Pval_neoblast_lower1 <- append(Wilcox_Pval_neoblast_lower1, w_nl1)
  
  w_u1l1 <- wdf$lower1[2]
  Wilcox_Pval_upper1_lower1 <- append(Wilcox_Pval_upper1_lower1, w_u1l1)
  
  w_u2l1 <- wdf$lower1[3]
  Wilcox_Pval_upper2_lower1 <- append(Wilcox_Pval_upper2_lower1, w_u2l1)
  
  w_u1n <- wdf$neoblast[2]
  Wilcox_Pval_upper1_neoblast <- append(Wilcox_Pval_upper1_neoblast, w_u1n)
  
  w_u2n <- wdf$neoblast[3]
  Wilcox_Pval_upper2_neoblast <- append(Wilcox_Pval_upper2_neoblast, w_u2n)
  
  w_u2u1 <- wdf$upper1[3]
  Wilcox_Pval_upper2_upper1 <- append(Wilcox_Pval_upper2_upper1, w_u2u1)
  
}


upper1 <- as.numeric(upper1)
upper1 <- round(upper1, digits = 2)

upper2 <- as.numeric(upper2)
upper2 <- round(upper2, digits = 2)

lower1 <- as.numeric(lower1)
lower1 <- round(lower1, digits = 2)

neoblast<- as.numeric(neoblast)
neoblast <- round(neoblast, digits = 2)

FvalANOVA<- as.numeric(FvalANOVA)
FvalANOVA <- round(FvalANOVA, digits = 0)

PvalANOVA<- as.numeric(PvalANOVA)
PvalANOVA<- round(PvalANOVA, digits = 4)

Tukey_Pval_neoblast_lower1<- as.numeric(Tukey_Pval_neoblast_lower1)
Tukey_Pval_neoblast_lower1<- round(Tukey_Pval_neoblast_lower1, digits = 4)

Tukey_Pval_upper1_lower1<- as.numeric(Tukey_Pval_upper1_lower1)
Tukey_Pval_upper1_lower1<- round(Tukey_Pval_upper1_lower1, digits = 4)

Tukey_Pval_upper2_lower1<- as.numeric(Tukey_Pval_upper2_lower1)
Tukey_Pval_upper2_lower1<- round(Tukey_Pval_upper2_lower1, digits = 4)

Tukey_Pval_upper1_neoblast<- as.numeric(Tukey_Pval_upper1_neoblast)
Tukey_Pval_upper1_neoblast<- round(Tukey_Pval_upper1_neoblast, digits = 4)

Tukey_Pval_upper2_neoblast<- as.numeric(Tukey_Pval_upper2_neoblast)
Tukey_Pval_upper2_neoblast<- round(Tukey_Pval_upper2_neoblast, digits = 4)

Tukey_Pval_upper2_upper1<- as.numeric(Tukey_Pval_upper2_upper1)
Tukey_Pval_upper2_upper1<- round(Tukey_Pval_upper2_upper1, digits = 4)

KruskalWallis<- as.numeric(KruskalWallis)
KruskalWallis<- round(KruskalWallis, digits = 0)

PvalKruskal<- as.numeric(PvalKruskal)
PvalKruskal<- round(PvalKruskal, digits = 4)

Wilcox_Pval_neoblast_lower1<- as.numeric(Wilcox_Pval_neoblast_lower1)
Wilcox_Pval_neoblast_lower1<- round(Wilcox_Pval_neoblast_lower1, digits = 4)

Wilcox_Pval_upper1_lower1<- as.numeric(Wilcox_Pval_upper1_lower1)
Wilcox_Pval_upper1_lower1<- round(Wilcox_Pval_upper1_lower1, digits = 4)

Wilcox_Pval_upper2_lower1<- as.numeric(Wilcox_Pval_upper2_lower1)
Wilcox_Pval_upper2_lower1<- round(Wilcox_Pval_upper2_lower1, digits = 4)

Wilcox_Pval_upper1_neoblast<- as.numeric(Wilcox_Pval_upper1_neoblast)
Wilcox_Pval_upper1_neoblast<- round(Wilcox_Pval_upper1_neoblast, digits = 4)

Wilcox_Pval_upper2_neoblast<- as.numeric(Wilcox_Pval_upper2_neoblast)
Wilcox_Pval_upper2_neoblast<- round(Wilcox_Pval_upper2_neoblast, digits = 4)

Wilcox_Pval_upper2_upper1<- as.numeric(Wilcox_Pval_upper2_upper1)
Wilcox_Pval_upper2_upper1<- round(Wilcox_Pval_upper2_upper1, digits = 4)

s_df <- rbind(names,
              ids,
              upper1,
              upper2,
              lower1,
              neoblast,
              FvalANOVA,
              PvalANOVA,
              Tukey_Pval_neoblast_lower1,
              Tukey_Pval_upper1_lower1,
              Tukey_Pval_upper1_neoblast,
              Tukey_Pval_upper2_lower1,
              Tukey_Pval_upper2_neoblast,
              Tukey_Pval_upper2_upper1,
              KruskalWallis, 
              PvalKruskal, 
              Wilcox_Pval_neoblast_lower1, 
              Wilcox_Pval_upper1_lower1, 
              Wilcox_Pval_upper1_neoblast,
              Wilcox_Pval_upper2_lower1,
              Wilcox_Pval_upper2_neoblast,
              Wilcox_Pval_upper2_upper1)

s_df <- data.frame(s_df)
s_df<- rotate_df(s_df)

#write_xlsx(s_df, "TegvsNeo_stats_table.xlsx")

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
genesteg <- genesteg[,-c(1,4:10)]
#write_xlsx(genesteg,"genesteg.xlsx")

#scale_y_discrete labels from opposite way to plotting so label from an inverted dataframe
library(purrr)

#try making heatmap read raw count expression slot

genestegflip <- genesteg
genestegflip <- genestegflip %>% map_df(rev)

write_xlsx(genestegflip, "genestegflip.xlsx") 

genestegflipSD <- read.csv("genestegflipSD.csv")
colnames(genestegflipSD) <- c('name', 'ID')


DoHeatmap(object = seurat_teg, features = genesteg$gene_ID, group.by = "cluster", assay = 'integrated', slot = "scale.data") +
  scale_fill_gradientn(colours = c('blue','black', 'orange')) +
  scale_y_discrete(labels = genestegflipSD$name)


####
#anova to compare expression of genes betwwen clusters
dat2 <- seurat_teg@assays[["integrated"]]@scale.data
dat2 <- data.frame(dat2)


#pick genes fro comparison 

num <- 1:66

fileConn<-file("original_cluster_stats_teg.txt")
sink(fileConn, append=TRUE)

for (i in num){
  
  n <- genesteg$protein[i]
  id <- genesteg$gene_ID[i]
  
  print(paste('Stats for gene', n, id))
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  tegcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(tegcountdf) <- c('ID', 'cell', 'expression')
  
  tegdata <- merge(tegcountdf, cluster_df, by = 'cell')
  
  #Mean expression for each group:
  print('Mean expression:')
  print(tegdata %>% 
          group_by(cluster) %>% 
          summarise(mean(expression)))  
  
  #make Graph
  my_comparisons = list(c("upper1", "lower1"), c("upper2","upper1"), c("upper2","lower1"))
  assign(paste("Plot", i, sep = ''),ggplot(data = tegdata, aes(x = cluster, y = expression)) +
           geom_violin() +
           labs(title = paste(n,"\n", id))+
           theme(plot.title = element_text(size=8))+
           geom_jitter(width=0.15, alpha=0.5) +
           stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = 'p.signif', paired = FALSE))


  #Anova to compare expression across the three clusters:
  mod <- aov(expression ~ cluster, data = tegdata)
  #View anova:
  print('Parametric ANOVA and TukeyHSD results:')
  print(summary(mod))
  print(TukeyHSD(mod))
  
  #Check normality:
  assign(paste("hist", i, sep = ''), as.ggplot(~qqPlot(mod$residuals)))
  #Check equal varience:
  assign(paste("box", i, sep = ''), as.ggplot(~boxplot(expression ~ cluster, data = tegdata)))
  
  
  print("Non-parametric - Kruskal Wallis and unpaired Wilcox test:")
  #If varience isnt equal: Welch's anova 
  print(kruskal.test(expression ~ cluster, data = tegdata))
  print(pairwise.wilcox.test(tegdata$expression, tegdata$cluster,
                             p.adjust.method = "none", paired = FALSE))
  
}
sink()
close(fileConn)

statplot <- ggarrange(Plot1,hist1,box1,    Plot2,hist2,box2,    Plot3,hist3,box3,    Plot4,hist4,box4,    Plot5,hist5,box5,    Plot6,hist6,box6,  Plot7,hist7,box7,    Plot8,hist8,box8,    Plot9,hist9,box9,    Plot10,hist10,box10, Plot11,hist11,box11, Plot12,hist12,box12, Plot13,hist13,box13, Plot14,hist14,box14, Plot15,hist15,box15, Plot16,hist16,box16, Plot17,hist17,box17, Plot18,hist18,box18, Plot19,hist19,box19, Plot20,hist20,box20, Plot21,hist21,box21, Plot22,hist22,box22, Plot23,hist23,box23, Plot24,hist24,box24, Plot25,hist25,box25, Plot26,hist26,box26, Plot27,hist27,box27, Plot28,hist28,box28, Plot29,hist29,box29, Plot30,hist30,box30, Plot31,hist31,box31, Plot32,hist32,box32, Plot33,hist33,box33, Plot34,hist34,box34, Plot35,hist35,box35, Plot36,hist36,box36, Plot37,hist37,box37, Plot38,hist38,box38, Plot39,hist39,box39, Plot40,hist40,box40, Plot41,hist41,box41, Plot42,hist42,box42, Plot43,hist43,box43, Plot44,hist44,box44, Plot45,hist45,box45, Plot46,hist46,box46, Plot47,hist47,box47, Plot48,hist48,box48, Plot49,hist49,box49, Plot50,hist50,box50, Plot51,hist51,box51, Plot52,hist52,box52, Plot53,hist53,box53, Plot54,hist54,box54, Plot55,hist55,box55, Plot56,hist56,box56, Plot57,hist57,box57, Plot58,hist58,box58, Plot59,hist59,box59, Plot60,hist60,box60, Plot61,hist61,box61, Plot62,hist62,box62, Plot63,hist63,box63, Plot64,hist64,box64, Plot65,hist65,box65, Plot66,hist66,box66, ncol=3, nrow=4)
#grpah saving
#####
png("original_cluster_stats_teg_1.png", width = 500, height = 700)
plot(statplot$`1`)
dev.off()

png("original_cluster_stats_teg_2.png", width = 500, height = 700)
plot(statplot$`2`)
dev.off()

png("original_cluster_stats_teg_3.png", width = 500, height = 700)
plot(statplot$`3`)
dev.off()

png("original_cluster_stats_teg_4.png", width = 500, height = 700)
plot(statplot$`4`)
dev.off()

png("original_cluster_stats_teg_5.png", width = 500, height = 700)
plot(statplot$`5`)
dev.off()

png("original_cluster_stats_teg_6.png", width = 500, height = 700)
plot(statplot$`6`)
dev.off()

png("original_cluster_stats_teg_7.png", width = 500, height = 700)
plot(statplot$`7`)
dev.off()

png("original_cluster_stats_teg_8.png", width = 500, height = 700)
plot(statplot$`8`)
dev.off()

png("original_cluster_stats_teg_9.png", width = 500, height = 700)
plot(statplot$`9`)
dev.off()

png("original_cluster_stats_teg_10.png", width = 500, height = 700)
plot(statplot$`10`)
dev.off()

png("original_cluster_stats_teg_11.png", width = 500, height = 700)
plot(statplot$`11`)
dev.off()

png("original_cluster_stats_teg_12.png", width = 500, height = 700)
plot(statplot$`12`)
dev.off()

png("original_cluster_stats_teg_13.png", width = 500, height = 700)
plot(statplot$`13`)
dev.off()

png("original_cluster_stats_teg_14.png", width = 500, height = 700)
plot(statplot$`14`)
dev.off()

png("original_cluster_stats_teg_15.png", width = 500, height = 700)
plot(statplot$`15`)
dev.off()

png("original_cluster_stats_teg_16.png", width = 500, height = 700)
plot(statplot$`16`)
dev.off()

png("original_cluster_stats_teg_17.png", width = 500, height = 700)
plot(statplot$`17`)
dev.off()
#####
names <- list()
ids <- list()
lower1 <- list()
neoblast <- list()
upper1 <- list()
upper2 <- list()
FvalANOVA <- list()
PvalANOVA <- list()
Tukey_Pval_upper1_lower1 <- list()
Tukey_Pval_upper2_lower1 <- list()
Tukey_Pval_upper2_upper1 <- list()
KruskalWallis <- list()
PvalKruskal <- list()
Wilcox_Pval_upper1_lower1 <- list()
Wilcox_Pval_upper2_lower1 <- list()
Wilcox_Pval_upper2_upper1 <- list()

for (i in num){
  
  n <- genesteg$protein[i]
  names <- append(names, n)
  id <- genesteg$gene_ID[i]
  ids <- append(ids, id)
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  tegcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(tegcountdf) <- c('ID', 'cell', 'expression')
  
  tegdata <- merge(tegcountdf, cluster_df, by = 'cell')
  
  m <- tegdata %>% 
          group_by(cluster) %>% 
          summarise(mean(expression))
  
  l1 <- m$`mean(expression)`[1]
  lower1 <- append(lower1, l1)
  u1 <- m$`mean(expression)`[2]
  upper1 <- append(upper1, u1)
  u2 <- m$`mean(expression)`[3]
  upper2 <- append(upper2, u2)
  

  mod <- aov(expression ~ cluster, data = tegdata)
  
  sum <- summary(mod)
  f <- sum[[1]][["F value"]][1]
  FvalANOVA <- append(f, FvalANOVA)
  p <- sum[[1]][["Pr(>F)"]][1]
  PvalANOVA <- append(p, PvalANOVA)
  
  tukey <- TukeyHSD(mod)
  tukeydf <- data.frame(tukey$cluster)
  
  u1_l1 <- tukeydf$p.adj[1]
  Tukey_Pval_upper1_lower1<- append(Tukey_Pval_upper1_lower1, u1_l1)
  u2_l1 <- tukeydf$p.adj[2]
  Tukey_Pval_upper2_lower1 <- append(Tukey_Pval_upper2_lower1, u2_l1)
  u2_u1 <- tukeydf$p.adj[3]
  Tukey_Pval_upper2_upper1 <- append(Tukey_Pval_upper2_upper1, u2_u1)

    #If varience isnt equal: Welch's anova 
  k <- kruskal.test(expression ~ cluster, data = tegdata)
  chi2 <- k[["statistic"]][["Kruskal-Wallis chi-squared"]]
  KruskalWallis <- append(KruskalWallis, chi2)
  pkw <- k[["p.value"]]
  PvalKruskal <- append(PvalKruskal, pkw)
  
  w <- pairwise.wilcox.test(tegdata$expression, tegdata$cluster,
                             p.adjust.method = "none", paired = FALSE)
  wdf <- data.frame(w$p.value)
  w_u1l1 <- wdf$lower1[1]
  Wilcox_Pval_upper1_lower1 <- append(Wilcox_Pval_upper1_lower1, w_u1l1)
  w_u2l1 <- wdf$lower1[2]
  Wilcox_Pval_upper2_lower1 <- append(Wilcox_Pval_upper2_lower1, w_u2l1)
  w_u2u1 <- wdf$upper1[2]
  Wilcox_Pval_upper2_upper1 <- append(Wilcox_Pval_upper2_upper1, w_u2u1)
}

upper1 <- as.numeric(upper1)
upper1 <- round(upper1, digits = 2)

upper2 <- as.numeric(upper2)
upper2 <- round(upper2, digits = 2)

lower1 <- as.numeric(lower1)
lower1 <- round(lower1, digits = 2)

FvalANOVA<- as.numeric(FvalANOVA)
FvalANOVA <- round(FvalANOVA, digits = 0)

PvalANOVA<- as.numeric(PvalANOVA)
PvalANOVA<- round(PvalANOVA, digits = 4)

Tukey_Pval_upper1_lower1<- as.numeric(Tukey_Pval_upper1_lower1)
Tukey_Pval_upper1_lower1<- round(Tukey_Pval_upper1_lower1, digits = 4)

Tukey_Pval_upper2_lower1<- as.numeric(Tukey_Pval_upper2_lower1)
Tukey_Pval_upper2_lower1<- round(Tukey_Pval_upper2_lower1, digits = 4)

Tukey_Pval_upper2_upper1<- as.numeric(Tukey_Pval_upper2_upper1)
Tukey_Pval_upper2_upper1<- round(Tukey_Pval_upper2_upper1, digits = 4)

KruskalWallis<- as.numeric(KruskalWallis)
KruskalWallis<- round(KruskalWallis, digits = 0)

PvalKruskal<- as.numeric(PvalKruskal)
PvalKruskal<- round(PvalKruskal, digits = 4)

Wilcox_Pval_upper1_lower1<- as.numeric(Wilcox_Pval_upper1_lower1)
Wilcox_Pval_upper1_lower1<- round(Wilcox_Pval_upper1_lower1, digits = 4)

Wilcox_Pval_upper2_lower1<- as.numeric(Wilcox_Pval_upper2_lower1)
Wilcox_Pval_upper2_lower1<- round(Wilcox_Pval_upper2_lower1, digits = 4)

Wilcox_Pval_upper2_upper1<- as.numeric(Wilcox_Pval_upper2_upper1)
Wilcox_Pval_upper2_upper1<- round(Wilcox_Pval_upper2_upper1, digits = 4)

s_df2 <- rbind(names,
              ids,
              upper1,
              upper2,
              lower1,
              FvalANOVA,
              PvalANOVA,
              Tukey_Pval_upper1_lower1,
              Tukey_Pval_upper2_lower1,
              Tukey_Pval_upper2_upper1,
              KruskalWallis, 
              PvalKruskal, 
              Wilcox_Pval_upper1_lower1, 
              Wilcox_Pval_upper2_lower1,
              Wilcox_Pval_upper2_upper1)

s_df2 <- data.frame(s_df2)
s_df2<- rotate_df(s_df2)

#write_xlsx(s_df2, "Tegcluster_stats_table.xlsx")

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
DimPlot(seurat_teg, reduction = "umap", pt.size = 2)

FeaturePlot(seurat_teg, features = 'Smp-051920',  slot = 'scale.data') + 
  scale_colour_gradientn(colours = c("#ABD9E9", "#E0F3F8", "#FDAE61", "#F46D43", "#D73027"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), name = paste('Expression', '\n', 'Z-score')) +
  labs(title = paste("nanos2","\n", "Smp-051920") )+
  theme(plot.title = element_text(size=10))

#SAVE PLOTS
genes <- c('Smp-051920','Smp-105360','Smp-175590','Smp-045200','Smp-346900')
names <- c('nanos-2','notch','fgfra','tal','sm25')
stat <- cbind(genes, names)
stat <- data.frame(stat)
num <- 1:5
#550 550
for (i in num){
  png(sprintf("%s.png", paste('reclusterteg', i)), width = 550, height = 550)
  plot(FeaturePlot(seurat_teg, features = stat$genes[i],  slot = 'data', pt.size = 2) + 
         labs(title = paste(stat$names[i],"\n", stat$genes[i]) )+
         theme(plot.title = element_text(size=10)))
  dev.off()
} 

for (i in num){
  png(sprintf("%s.png", paste('reclustertegvln', i)), width = 850, height = 500)
  plot(VlnPlot(seurat_teg, features = stat$genes[i], pt.size = 0.2) +
         theme(legend.position = "none", plot.title = element_text(size=10))+
         labs(title = paste(stat$names[i],"\n", stat$genes[i])))
  dev.off()
}
#550 550


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


write_xlsx(genestegflip, "genestegflip.xlsx") 


genestegflipRCSD <- read.csv("genestegflipRCSD.csv")
colnames(genestegflipRCSD) <- c('name', 'ID')


#heatmap for tegument proteins on reclustered cells
DoHeatmap(object = seurat_teg, features = genesteg$gene_ID, group.by = "seurat_clusters", assay = 'RNA', slot = "scale.data") +
  scale_fill_gradientn(colours = c('blue', 'black', 'orange')) +
  scale_y_discrete(labels = genestegflipRCSD$name)

#stats on recluster


#pick genes fro comparison 

num <- 1:66

dat2 <- seurat_teg@assays[["RNA"]]@scale.data
dat2 <- data.frame(dat2)
teg <- seurat_teg@active.ident
tegdf <- data.frame(teg)
setDT(tegdf, keep.rownames = TRUE)[]
colnames(tegdf) <- c("cell", "cluster")
tegdf$cluster <- as.factor(tegdf$cluster)

fileConn<-file("recluster_stats_teg.txt")
sink(fileConn, append=TRUE)

for (i in num){
  
  n <- genesteg$protein[i]
  id <- genesteg$gene_ID[i]
  
  print(paste('Stats for gene', n, id))
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  tegcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(tegcountdf) <- c('ID', 'cell', 'expression')
  
  tegdata <- merge(tegcountdf, tegdf, by = 'cell')
  
  #Mean expression for each group:
  print('Mean expression:')
  print(tegdata %>% 
          group_by(cluster) %>% 
          summarise(mean(expression)))  
  
  assign(paste("Plot", i, sep = ''),ggplot(data = tegdata, aes(x = cluster, y = expression)) +
           geom_violin() +
           labs(title = paste(n,"\n", id))+
           theme(plot.title = element_text(size=8))+
           geom_jitter(width=0.15, alpha=0.5))

  #Anova to compare expression across the three clusters:
  mod <- aov(expression ~ cluster, data = tegdata)
  #View anova:
  print('Parametric ANOVA and TukeyHSD results:')
  print(summary(mod))
  print(TukeyHSD(mod))
  
  #Check normality:
  assign(paste("hist", i, sep = ''), as.ggplot(~qqPlot(mod$residuals)))
  #Check equal varience:
  assign(paste("box", i, sep = ''), as.ggplot(~boxplot(expression ~ cluster, data = tegdata)))
  
  
  print("Non-parametric - Kruskal Wallis and unpaired Wilcox test:")
  print(kruskal.test(expression ~ cluster, data = tegdata))
  print(pairwise.wilcox.test(tegdata$expression, tegdata$cluster,
                             p.adjust.method = "none", paired = FALSE))
  
}
sink()
close(fileConn)

statplot <- ggarrange(Plot1,hist1,box1,    Plot2,hist2,box2,    Plot3,hist3,box3,    Plot4,hist4,box4,    Plot5,hist5,box5,    Plot6,hist6,box6,  Plot7,hist7,box7,    Plot8,hist8,box8,    Plot9,hist9,box9,    Plot10,hist10,box10, Plot11,hist11,box11, Plot12,hist12,box12, Plot13,hist13,box13, Plot14,hist14,box14, Plot15,hist15,box15, Plot16,hist16,box16, Plot17,hist17,box17, Plot18,hist18,box18, Plot19,hist19,box19, Plot20,hist20,box20, Plot21,hist21,box21, Plot22,hist22,box22, Plot23,hist23,box23, Plot24,hist24,box24, Plot25,hist25,box25, Plot26,hist26,box26, Plot27,hist27,box27, Plot28,hist28,box28, Plot29,hist29,box29, Plot30,hist30,box30, Plot31,hist31,box31, Plot32,hist32,box32, Plot33,hist33,box33, Plot34,hist34,box34, Plot35,hist35,box35, Plot36,hist36,box36, Plot37,hist37,box37, Plot38,hist38,box38, Plot39,hist39,box39, Plot40,hist40,box40, Plot41,hist41,box41, Plot42,hist42,box42, Plot43,hist43,box43, Plot44,hist44,box44, Plot45,hist45,box45, Plot46,hist46,box46, Plot47,hist47,box47, Plot48,hist48,box48, Plot49,hist49,box49, Plot50,hist50,box50, Plot51,hist51,box51, Plot52,hist52,box52, Plot53,hist53,box53, Plot54,hist54,box54, Plot55,hist55,box55, Plot56,hist56,box56, Plot57,hist57,box57, Plot58,hist58,box58, Plot59,hist59,box59, Plot60,hist60,box60, Plot61,hist61,box61, Plot62,hist62,box62, Plot63,hist63,box63, Plot64,hist64,box64, Plot65,hist65,box65, Plot66,hist66,box66, ncol=3, nrow=3)
#graph saving
#####
png("recluster_statsplot_teg_1.png", width = 470, height = 800)
plot(statplot$`1`)
dev.off()

png("recluster_statsplot_teg_2.png", width = 470, height = 800)
plot(statplot$`2`)
dev.off()


png("recluster_statsplot_teg_3.png", width = 470, height = 800)
plot(statplot$`3`)
dev.off()

png("recluster_statsplot_teg_4.png", width = 470, height = 800)
plot(statplot$`4`)
dev.off()

png("recluster_statsplot_teg_5.png", width = 470, height = 800)
plot(statplot$`5`)
dev.off()

png("recluster_statsplot_teg_6.png", width = 470, height = 800)
plot(statplot$`6`)
dev.off()

png("recluster_statsplot_teg_7.png", width = 470, height = 800)
plot(statplot$`7`)
dev.off()

png("recluster_statsplot_teg_8.png", width = 470, height = 800)
plot(statplot$`8`)
dev.off()

png("recluster_statsplot_teg_9.png", width = 470, height = 800)
plot(statplot$`9`)
dev.off()

png("recluster_statsplot_teg_10.png", width = 470, height = 800)
plot(statplot$`10`)
dev.off()

png("recluster_statsplot_teg_11.png", width = 470, height = 800)
plot(statplot$`11`)
dev.off()

png("recluster_statsplot_teg_12.png", width = 470, height = 800)
plot(statplot$`12`)
dev.off()

png("recluster_statsplot_teg_13.png", width = 470, height = 800)
plot(statplot$`13`)
dev.off()

png("recluster_statsplot_teg_14.png", width = 470, height = 800)
plot(statplot$`14`)
dev.off()

png("recluster_statsplot_teg_15.png", width = 470, height = 800)
plot(statplot$`15`)
dev.off()

png("recluster_statsplot_teg_16.png", width = 470, height = 800)
plot(statplot$`16`)
dev.off()

png("recluster_statsplot_teg_17.png", width = 470, height = 800)
plot(statplot$`17`)
dev.off()

png("recluster_statsplot_teg_18.png", width = 470, height = 800)
plot(statplot$`18`)
dev.off()

png("recluster_statsplot_teg_19.png", width = 470, height = 800)
plot(statplot$`19`)
dev.off()

png("recluster_statsplot_teg_20.png", width = 470, height = 800)
plot(statplot$`20`)
dev.off()

png("recluster_statsplot_teg_21.png", width = 470, height = 800)
plot(statplot$`21`)
dev.off()

png("recluster_statsplot_teg_22.png", width = 470, height = 800)
plot(statplot$`22`)
dev.off()
## stats to table
names <- list()
ids <- list()
FvalANOVA <- list()
PvalANOVA <- list()
KruskalWallis <- list()
PvalKruskal <- list()
o0 <- list()
i1 <- list()
ii2 <- list()
iii3 <- list()
iv4 <- list()
v5 <- list()
vi6 <- list()
vii7 <- list()

for (i in num){
  
  n <- genesteg$protein[i]
  names <- append(names, n)
  id <- genesteg$gene_ID[i]
  ids <- append(ids, id)
  
  newdat2 <- dat2[id,]
  cn2 <- colnames(newdat2)
  setDT(newdat2, keep.rownames = TRUE)[]
  tegcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  colnames(tegcountdf) <- c('ID', 'cell', 'expression')
  tegdata <- merge(tegcountdf, tegdf, by = 'cell')
  
  #Mean expression for each group:
   m <- tegdata %>% 
           group_by(cluster) %>% 
           summarise(mean(expression))
   

   o <- m$`mean(expression)`[1]
   o0 <- append(o0, o)
   i <- m$`mean(expression)`[2]
   i1 <- append(i1, i)
   ii <- m$`mean(expression)`[3]
   ii2 <- append(ii2, ii)
   iii <- m$`mean(expression)`[4]
   iii3 <- append(iii3, iii)
   iv <- m$`mean(expression)`[5]
   iv4 <- append(iv4, iv)
   v <- m$`mean(expression)`[6]
   v5 <- append(v5, v)
   vi <- m$`mean(expression)`[7]
   vi6 <- append(vi6, vi)
   vii <- m$`mean(expression)`[8]
   vii7 <- append(vii7, vii)
   
   
   
  # 

  #Anova to compare expression across the three clusters:
  mod <- aov(expression ~ cluster, data = tegdata)
  #View anova:
  sum <- summary(mod)
  f <- sum[[1]][["F value"]][1]
  FvalANOVA <- append(f, FvalANOVA)
  p <- sum[[1]][["Pr(>F)"]][1]
  PvalANOVA <- append(p, PvalANOVA)
  
  #tukey <- TukeyHSD(mod)
  
  #Check normality:
  k <- kruskal.test(expression ~ cluster, data = tegdata)
  chi2 <- k[["statistic"]][["Kruskal-Wallis chi-squared"]]
  KruskalWallis <- append(KruskalWallis, chi2)
  pkw <- k[["p.value"]]
  PvalKruskal <- append(PvalKruskal, pkw)
  
  #w <- pairwise.wilcox.test(tegdata$expression, tegdata$cluster,
  #                           p.adjust.method = "none", paired = FALSE)
  
  
}

FvalANOVA<- as.numeric(FvalANOVA)
FvalANOVA <- round(FvalANOVA, digits = 0)

PvalANOVA<- as.numeric(PvalANOVA)
PvalANOVA<- round(PvalANOVA, digits = 4)

KruskalWallis<- as.numeric(KruskalWallis)
KruskalWallis<- round(KruskalWallis, digits = 0)

PvalKruskal<- as.numeric(PvalKruskal)
PvalKruskal<- round(PvalKruskal, digits = 4)



s_df3 <- rbind(names,
               ids,
               FvalANOVA,
               PvalANOVA,
               KruskalWallis, 
               PvalKruskal,
               o0,
               i1,
               ii2,
               iii3,
               iv4,
               v5,
               vi6,
               vii7)

s_df3 <- data.frame(s_df3)
s_df3<- rotate_df(s_df3)

#write_xlsx(s_df3, "TegRecluster_stats_table.xlsx")
#####

#finds markers that disinguish the sub clusters from each other
allmarkersT <- FindAllMarkers(seurat_teg, assay = 'RNA', test.use =  'roc', only.pos = TRUE, return.thresh = 0)
markersT <- filter(allmarkersT, myAUC >= 0.8)

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






