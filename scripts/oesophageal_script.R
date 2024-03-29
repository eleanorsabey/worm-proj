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
library(sjmisc)

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


#Plots 
#####
#Map of whole Seurat object:
DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + theme(legend.position = "none") 
#Oesophageal gland cells highlighted amongst all cells:
DimPlot(object = seurat_object, reduction = 'umap', cells.highlight = c) + theme(legend.position = "none") 
#Plot only oesophageal gland cells:
DimPlot(object = seurat_object, cells = c, reduction = 'umap')
#####

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


#pick genes fro comparison 
genes <- c('Smp-051920','Smp-105360','Smp-175590','Smp-125320','Smp-010550', 'Smp-172180', 'Smp-085840', 'Smp-124000')
names <- c('nanos-2','notch','fgfra','meg-9','meg-15', 'meg-8.2', 'meg-4.2','meg-14')
stat <- cbind(genes, names)
stat <- data.frame(stat)
num <- 1:8


dat2 <- seurat_og_neo@assays[["integrated"]]@scale.data
dat2 <- data.frame(dat2)


fileConn<-file("og_vs_neo_stats.txt")
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
  assign(paste("hist", i, sep = ''), as.ggplot(~qqPlot(mod$residuals)))
  
  #Check equal varience:
  assign(paste("box", i, sep = ''), as.ggplot(~boxplot(expression ~ cluster, data = ogandneo)))
  

  print("Non-parametric - Kruskal Wallis and unpaired Wilcox test:")
  #If varience isnt equal: Welch's anova 
  print(kruskal.test(expression ~ cluster, data = ogandneo))
  print(pairwise.wilcox.test(ogandneo$expression, ogandneo$cluster,
                       p.adjust.method = "none", paired = FALSE))
  
}
sink()
close(fileConn)
statplot <- ggarrange(Plot1, hist1, box1, Plot2, hist2, box2, Plot3, hist3, box3, Plot4,  hist4, box4, ncol=3, nrow=4)
statplot2 <- ggarrange(Plot5, hist5, box5, Plot6,  hist6, box6, Plot7, hist7, box7, Plot8, hist8, box8, ncol=3, nrow=4)
plot(statplot)
png("og_vs_neo_stats_1.png", width = 470, height = 800)
plot(statplot)
dev.off()
png("og_vs_neo_stats_2.png", width = 470, height = 800)
plot(statplot2)
dev.off()
  
## stats loop that writes a table
names <- list()
ids <- list()
lower <- list()
neoblast <- list()
upper <- list()
FvalANOVA <- list()
PvalANOVA <- list()
Tukey_Pval_neoblast_lower <- list()
Tukey_Pval_upper_lower <- list()
Tukey_Pval_upper_neoblast <- list()
KruskalWallis <- list()
PvalKruskal <- list()
Wilcox_Pval_neoblast_lower <- list()
Wilcox_Pval_upper_lower <- list()
Wilcox_Pval_upper_neoblast <- list()

for (i in num){
  
  n <- stat$names[i]
  names <- append(names, n)
  id <- stat$genes[i]
  ids <- append(ids, id)

  # 
   newdat2 <- dat2[id,]
  # 
   cn2 <- colnames(newdat2)
  # 
   setDT(newdat2, keep.rownames = TRUE)[]
  # 
   ognbcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  # 
   colnames(ognbcountdf) <- c('ID', 'cell', 'expression')
  # 
   ogandneo <- merge(ognbcountdf, neo_og_cells, by = 'cell')
  # 

   m <- ogandneo %>% 
          group_by(cluster) %>% 
          summarise(mean(expression))
   
  l <- m$`mean(expression)`[1]
  lower <- append(lower, l)
  n <- m$`mean(expression)`[2]
  neoblast <- append(neoblast, n)
  u <- m$`mean(expression)`[3]
  upper <- append(upper,u)

   mod <- aov(expression ~ cluster, data = ogandneo)

   sum <- summary(mod)
   f <- sum[[1]][["F value"]][1]
   FvalANOVA <- append(f, FvalANOVA)
   p <- sum[[1]][["Pr(>F)"]][1]
   PvalANOVA <- append(p, PvalANOVA)
   
  tukey <- TukeyHSD(mod)
  tukeydf <- data.frame(tukey$cluster)
  n_l <- tukeydf$p.adj[1]
  Tukey_Pval_neoblast_lower<- append(Tukey_Pval_neoblast_lower, n_l)
  u_l <- tukeydf$p.adj[2]
  Tukey_Pval_upper_lower <- append(Tukey_Pval_upper_lower, u_l)
  u_n <- tukeydf$p.adj[3]
  Tukey_Pval_upper_neoblast <- append(Tukey_Pval_upper_neoblast, u_n)
 
    k <- kruskal.test(expression ~ cluster, data = ogandneo)
    chi2 <- k[["statistic"]][["Kruskal-Wallis chi-squared"]]
    KruskalWallis <- append(KruskalWallis, chi2)
    pkw <- k[["p.value"]]
    PvalKruskal <- append(PvalKruskal, pkw)
    
   w <- pairwise.wilcox.test(ogandneo$expression, ogandneo$cluster,
                              p.adjust.method = "none", paired = FALSE)
  
  
   wdf <- data.frame(w$p.value)
   w_nl <- wdf$lower[1]
   Wilcox_Pval_neoblast_lower <- append(Wilcox_Pval_neoblast_lower, w_nl)
   w_ul <- wdf$lower[2]
   Wilcox_Pval_upper_lower <- append(Wilcox_Pval_upper_lower, w_ul)
   w_un <- wdf$neoblast[2]
   Wilcox_Pval_upper_neoblast <- append(Wilcox_Pval_upper_neoblast, w_un)
   
   up <- filter(ogandneo, ogandneo$cluster == "upper")
   lo <- filter(ogandneo, ogandneo$cluster == "lower")
   neo <- filter(ogandneo, ogandneo$cluster == "neoblast")
   
   wuplo <- wilcox.test(up$expression, lo$expression)
   wilcox_Wval_upper_lower <- append(wilcox_Wval_upper_lower, wuplo[["statistic"]])
   wupneo <- wilcox.test(up$expression, neo$expression)
   wloneo <- wilcox.test(lo$expression, neo$expression)
  
}


upper <- as.numeric(upper)
upper <- round(upper, digits = 2)

lower <- as.numeric(lower)
lower <- round(lower, digits = 2)

neoblast<- as.numeric(neoblast)
neoblast <- round(neoblast, digits = 2)

FvalANOVA<- as.numeric(FvalANOVA)
FvalANOVA <- round(FvalANOVA, digits = 0)

PvalANOVA<- as.numeric(PvalANOVA)
PvalANOVA<- round(PvalANOVA, digits = 4)

Tukey_Pval_neoblast_lower<- as.numeric(Tukey_Pval_neoblast_lower)
Tukey_Pval_neoblast_lower<- round(Tukey_Pval_neoblast_lower, digits = 4)

Tukey_Pval_upper_lower<- as.numeric(Tukey_Pval_upper_lower)
Tukey_Pval_upper_lower<- round(Tukey_Pval_upper_lower, digits = 4)

Tukey_Pval_upper_neoblast<- as.numeric(Tukey_Pval_upper_neoblast)
Tukey_Pval_upper_neoblast<- round(Tukey_Pval_upper_neoblast, digits = 4)

KruskalWallis<- as.numeric(KruskalWallis)
KruskalWallis<- round(KruskalWallis, digits = 0)

PvalKruskal<- as.numeric(PvalKruskal)
PvalKruskal<- round(PvalKruskal, digits = 4)

Wilcox_Pval_neoblast_lower<- as.numeric(Wilcox_Pval_neoblast_lower)
Wilcox_Pval_neoblast_lower<- round(Wilcox_Pval_neoblast_lower, digits = 4)

Wilcox_Pval_upper_lower<- as.numeric(Wilcox_Pval_upper_lower)
Wilcox_Pval_upper_lower<- round(Wilcox_Pval_upper_lower, digits = 4)

Wilcox_Pval_upper_neoblast<- as.numeric(Wilcox_Pval_upper_neoblast)
Wilcox_Pval_upper_neoblast<- round(Wilcox_Pval_upper_neoblast, digits = 4)

s_df <- rbind(names,
              ids,
              upper,
              lower,
              neoblast,
              FvalANOVA,
              PvalANOVA,
              Tukey_Pval_neoblast_lower,
              Tukey_Pval_upper_lower,
              Tukey_Pval_upper_neoblast,
              KruskalWallis, 
              PvalKruskal, 
              Wilcox_Pval_neoblast_lower, 
              Wilcox_Pval_upper_lower, 
              Wilcox_Pval_upper_neoblast)

s_df <- data.frame(s_df)
s_df<- rotate_df(s_df)

#write.csv(s_df, "OgvsNeo_stats_table.csv")

##

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

#write_xlsx(eso_proteinsflip, "eso_proteinsflip.xlsx") 

eso_proteinsflipSD <- read.csv("eso_proteinsflipSD.csv")
colnames(eso_proteinsflipSD) <- c('name', 'ID')
#Heatmap using scale.data
DoHeatmap(object = seurat_sub, group.by= 'cluster', features = eso_proteins$ID, draw.lines = TRUE)+ 
    scale_fill_gradientn(colours = c('blue', 'black', 'orange')) +
  scale_y_discrete(labels = eso_proteinsflipSD$name)


#It is decided that the upper cluster resembles oesophageal cells more than neoblast so will not be removed from the data before reclustering. 
#The outlier cell is also not removed becase heatmap expression fits with lower cluster. 


#t test to compare expression of genes betwwen upper and lower cluster



#pick genes fro comparison 

num <- 1:71

dat2 <- seurat_sub@assays[["integrated"]]@scale.data
dat2 <- data.frame(dat2)

fileConn<-file("original_cluster_stats_OG.txt")
sink(fileConn, append=TRUE)


for (i in num){
  
  n <- eso_proteins$name[i]
  id <- eso_proteins$ID[i]
  
  print(paste('Stats for gene', n, id))
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  ogcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(ogcountdf) <- c('ID', 'cell', 'expression')
  
  ogdata <- merge(ogcountdf, cluster_df, by = 'cell')
  
  #Mean expression for each group:
  print('Mean expression:')
  print(ogdata %>% 
          group_by(cluster) %>% 
          summarise(mean(expression)))  
  
  #make Graph
  my_comparisons = list(c("lower","upper"))
  assign(paste("Plot", i, sep = ''),ggplot(data = ogdata, aes(x = cluster, y = expression)) +
           geom_violin() +
           labs(title = paste(n,"\n", id))+
           theme(plot.title = element_text(size=8))+
           geom_jitter(width=0.15, alpha=0.5) +
           stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = 'p.signif', paired = FALSE))
  
  #check normaility 
  
  ogupper <- filter(ogdata, cluster == "upper")
  oglower <- filter(ogdata, cluster == "lower")
  assign(paste("histu", i, sep = ''), as.ggplot(~qqPlot(ogupper$expression)))
  assign(paste("histl", i, sep = ''), as.ggplot(~qqPlot(oglower$expression)))
  
  
  #check equal varience 
  print('Check equal varience:')
  print(var.test(expression ~ cluster, data = ogdata))
  
  #t test to compare expression across the clusters:
  print('Parametric Independent t-test:')
  print(t.test(data = ogdata,
         expression ~ cluster,
         var.equal = T))
  
  #non parametric alternative:
  print('Non-parametric unpaired two-sample Wilcoxon test:')
  print(wilcox.test(expression ~ cluster, data = ogdata))
  
}
sink()
close(fileConn)


statplot <- ggarrange(Plot1,histu1,histl1, Plot2,histu2,histl2, Plot3,histu3,histl3, Plot4,histu4,histl4, Plot5,histu5,histl5,  Plot6,histu6,histl6, Plot7,histu7,histl7, Plot8,histu8,histl8, Plot9,histu9,histl9, Plot10,histu10,histl10, Plot11,histu11,histl11, Plot12,histu12,histl12, Plot13,histu13,histl13, Plot14,histu14,histl14, Plot15,histu15,histl15,Plot16,histu16,histl16, Plot17,histu17,histl17, Plot18,histu18,histl18, Plot19,histu19,histl19, Plot20,histu20,histl20,Plot21,histu21,histl21, Plot22,histu22,histl22, Plot23,histu23,histl23, Plot24,histu24,histl24, Plot25,histu25,histl25,Plot26,histu26,histl26, Plot27,histu27,histl27, Plot28,histu28,histl28, Plot29,histu29,histl29, Plot30,histu30,histl30, Plot31,histu31,histl31, Plot32,histu32,histl32, Plot33,histu33,histl33, Plot34,histu34,histl34, Plot35,histu35,histl35, Plot36,histu36,histl36, Plot37,histu37,histl37, Plot38,histu38,histl38, Plot39,histu39,histl39, Plot40,histu40,histl40, Plot41,histu41,histl41, Plot42,histu42,histl42, Plot43,histu43,histl43, Plot44,histu44,histl44, Plot45,histu45,histl45, Plot46,histu46,histl46, Plot47,histu47,histl47, Plot48,histu48,histl48, Plot49,histu49,histl49, Plot50,histu50,histl50, Plot51,histu51,histl51, Plot52,histu52,histl52, Plot53,histu53,histl53, Plot54,histu54,histl54, Plot55,histu55,histl55, Plot56,histu56,histl56, Plot57,histu57,histl57, Plot58,histu58,histl58, Plot59,histu59,histl59, Plot60,histu60,histl60, Plot61,histu61,histl61, Plot62,histu62,histl62, Plot63,histu63,histl63, Plot64,histu64,histl64, Plot65,histu65,histl65, Plot66,histu66,histl66, Plot67,histu67,histl67, Plot68,histu68,histl68, Plot69,histu69,histl69, Plot70,histu70,histl70, Plot71,histu71,histl71,ncol=3, nrow=4)
#graph saving:
######
png("original_cluster_stats_OG_1.png", width = 470, height = 800)
plot(statplot$`1`)
dev.off()

png("original_cluster_stats_OG_2.png", width = 470, height = 800)
plot(statplot$`2`)
dev.off()

png("original_cluster_stats_OG_3.png", width = 470, height = 800)
plot(statplot$`3`)
dev.off()

png("original_cluster_stats_OG_4.png", width = 470, height = 800)
plot(statplot$`4`)
dev.off()

png("original_cluster_stats_OG_5.png", width = 470, height = 800)
plot(statplot$`5`)
dev.off()

png("original_cluster_stats_OG_6.png", width = 470, height = 800)
plot(statplot$`6`)
dev.off()

png("original_cluster_stats_OG_7.png", width = 470, height = 800)
plot(statplot$`7`)
dev.off()

png("original_cluster_stats_OG_8.png", width = 470, height = 800)
plot(statplot$`8`)
dev.off()

png("original_cluster_stats_OG_9.png", width = 470, height = 800)
plot(statplot$`9`)
dev.off()

png("original_cluster_stats_OG_10.png", width = 470, height = 800)
plot(statplot$`10`)
dev.off()

png("original_cluster_stats_OG_11.png", width = 470, height = 800)
plot(statplot$`11`)
dev.off()

png("original_cluster_stats_OG_12.png", width = 470, height = 800)
plot(statplot$`12`)
dev.off()

png("original_cluster_stats_OG_13.png", width = 470, height = 800)
plot(statplot$`13`)
dev.off()

png("original_cluster_stats_OG_14.png", width = 470, height = 800)
plot(statplot$`14`)
dev.off()

png("original_cluster_stats_OG_15.png", width = 470, height = 800)
plot(statplot$`15`)
dev.off()

png("original_cluster_stats_OG_16.png", width = 470, height = 800)
plot(statplot$`16`)
dev.off()

png("original_cluster_stats_OG_17.png", width = 470, height = 800)
plot(statplot$`17`)
dev.off()

png("original_cluster_stats_OG_18.png", width = 470, height = 800)
plot(statplot$`18`)
dev.off()

#####
#write stats into table
names <- list()
ids <- list()
lower <- list()
upper <- list()
Tval <- list()
PvalTtest <- list()
Wilcoxon_Pval <- list()

for (i in num){
  
  n <- eso_proteins$name[i]
  names <- append(names, n)
  id <- eso_proteins$ID[i]
  ids <- append(ids, id)
  #print(paste('Stats for gene', n, id))
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  ogcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(ogcountdf) <- c('ID', 'cell', 'expression')
  
  ogdata <- merge(ogcountdf, cluster_df, by = 'cell')
  
  #Mean expression for each group:
  #print('Mean expression:')
     m <- ogdata %>% 
            group_by(cluster) %>% 
            summarise(mean(expression))  
     
     l <- m$`mean(expression)`[1]
     lower <- append(lower, l)
     u <- m$`mean(expression)`[2]
     upper <- append(upper,u)
     
  ttest <- t.test(data = ogdata,
               expression ~ cluster,
               var.equal = T)
  tv <- ttest[["statistic"]][["t"]]
  Tval <- append(tv, Tval)
  pv <- ttest[["p.value"]]
  PvalTtest <- append(pv, PvalTtest)
  
  #non parametric alternative:
  #print('Non-parametric unpaired two-sample Wilcoxon test:')
  w <- wilcox.test(expression ~ cluster, data = ogdata)
  wdf <- data.frame(w$p.value)
  w_ul <- wdf$w.p.value[1]
  Wilcoxon_Pval <- append(Wilcoxon_Pval, w_ul)
}

upper <- as.numeric(upper)
upper <- round(upper, digits = 2)

lower <- as.numeric(lower)
lower <- round(lower, digits = 2)

Tval<- as.numeric(Tval)
Tval <- round(Tval, digits = 0)

PvalTtest<- as.numeric(PvalTtest)
PvalTtest<- round(PvalTtest, digits = 4)

Wilcoxon_Pval<- as.numeric(Wilcoxon_Pval)
Wilcoxon_Pval<- round(Wilcoxon_Pval, digits = 4)

s_df2 <- rbind(names,
              ids,
              upper,
              lower,
              Tval,
              PvalTtest,
              Wilcoxon_Pval
              )

s_df2 <- data.frame(s_df2)
s_df2<- rotate_df(s_df2)

#write_xlsx(s_df2, "OGcluster_stat_table.xlsx")

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

#write_xlsx(eso_proteinsflip, "eso_proteinsflip.xlsx") 

eso_proteinsflipRCSD <- read.csv("eso_proteinsflipRCSD.csv")
colnames(eso_proteinsflipRCSD) <- c('name', 'ID')

DoHeatmap(object = seurat_sub, features = eso_proteins$ID, group.by = "seurat_clusters", assay = 'RNA', slot = "scale.data", group.colors = c('#7AD169', '#928DE0', '#FAA637')) +
  scale_fill_gradientn(colours = c('blue', 'black', 'orange'), breaks=c(-3,-2,-1,0,1,2,3),labels=c(-3,-2,-1,0,1,2,3),limits=c(-2.5,2.5)) +
  scale_y_discrete(labels = eso_proteinsflipRCSD$name) 


#Compare expression of specific proteins across new clusters:


#pick genes fro comparison 

dat2 <- seurat_sub@assays[["RNA"]]@scale.data
dat2 <- data.frame(dat2)
og <- seurat_sub@active.ident
ogdf <- data.frame(og)
setDT(ogdf, keep.rownames = TRUE)[]
colnames(ogdf) <- c("cell", "cluster")
ogdf$cluster <- as.factor(ogdf$cluster)

fileConn<-file("recluster_stats_OG.txt")
sink(fileConn, append=TRUE)

num <- 1:71
for (i in num){
  
  n <- eso_proteins$name[i]
  id <- eso_proteins$ID[i]
  
  print(paste('Stats for gene', n, id))
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  ogcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(ogcountdf) <- c('ID', 'cell', 'expression')
  
  ogdata <- merge(ogcountdf, ogdf, by = 'cell')
  
  #Mean expression for each group:
  print('Mean expression:')
  print(ogdata %>% 
          group_by(cluster) %>% 
          summarise(mean(expression)))  
  
  my_comparisons = list(c('0', '1'), c('1', '2'), c('0', '2'))
  assign(paste("Plot", i, sep = ''),ggplot(data = ogdata, aes(x = cluster, y = expression)) +
           geom_violin() +
           labs(title = paste(n,"\n", id))+
           theme(plot.title = element_text(size=8))+
           geom_jitter(width=0.15, alpha=0.5) +
           stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = 'p.signif', paired = FALSE))
  
  #Anova to compare expression across the three clusters:
  mod <- aov(expression ~ cluster, data = ogdata)
  #View anova:
  print('Parametric ANOVA and TukeyHSD results:')
  print(summary(mod))
  print(TukeyHSD(mod))
  
  #Check normality:
  assign(paste("hist", i, sep = ''), as.ggplot(~qqPlot(mod$residuals)))
  #Check equal varience:
  assign(paste("box", i, sep = ''), as.ggplot(~boxplot(expression ~ cluster, data = ogdata)))
  
  
  print("Non-parametric - Kruskal Wallis and unpaired Wilcox test:")
  #If varience isnt equal: Welch's anova 
  print(kruskal.test(expression ~ cluster, data = ogdata))
  print(pairwise.wilcox.test(ogdata$expression, ogdata$cluster,
                             p.adjust.method = "none", paired = FALSE))
}

sink()
close(fileConn)


print(noquote(paste("Plot", 1:71, ",", "hist", 1:71, ",", "box", 1:71, ",", sep = '')))
statplot <- ggarrange(Plot1,hist1,box1,    Plot2,hist2,box2,    Plot3,hist3,box3,    Plot4,hist4,box4,    Plot5,hist5,box5,    Plot6,hist6,box6,  Plot7,hist7,box7,    Plot8,hist8,box8,    Plot9,hist9,box9,    Plot10,hist10,box10, Plot11,hist11,box11, Plot12,hist12,box12, Plot13,hist13,box13, Plot14,hist14,box14, Plot15,hist15,box15, Plot16,hist16,box16, Plot17,hist17,box17, Plot18,hist18,box18, Plot19,hist19,box19, Plot20,hist20,box20, Plot21,hist21,box21, Plot22,hist22,box22, Plot23,hist23,box23, Plot24,hist24,box24, Plot25,hist25,box25, Plot26,hist26,box26, Plot27,hist27,box27, Plot28,hist28,box28, Plot29,hist29,box29, Plot30,hist30,box30, Plot31,hist31,box31, Plot32,hist32,box32, Plot33,hist33,box33, Plot34,hist34,box34, Plot35,hist35,box35, Plot36,hist36,box36, Plot37,hist37,box37, Plot38,hist38,box38, Plot39,hist39,box39, Plot40,hist40,box40, Plot41,hist41,box41, Plot42,hist42,box42, Plot43,hist43,box43, Plot44,hist44,box44, Plot45,hist45,box45, Plot46,hist46,box46, Plot47,hist47,box47, Plot48,hist48,box48, Plot49,hist49,box49, Plot50,hist50,box50, Plot51,hist51,box51, Plot52,hist52,box52, Plot53,hist53,box53, Plot54,hist54,box54, Plot55,hist55,box55, Plot56,hist56,box56, Plot57,hist57,box57, Plot58,hist58,box58, Plot59,hist59,box59, Plot60,hist60,box60, Plot61,hist61,box61, Plot62,hist62,box62, Plot63,hist63,box63, Plot64,hist64,box64, Plot65,hist65,box65, Plot66,hist66,box66, Plot67,hist67,box67, Plot68,hist68,box68, Plot69,hist69,box69, Plot70,hist70,box70, Plot71,hist71,box71, ncol=3, nrow=4)

#graph saving
#####
png("recluster_statsplot_OG_1.png", width = 470, height = 800)
plot(statplot$`1`)
dev.off()

png("recluster_statsplot_OG_2.png", width = 470, height = 800)
plot(statplot$`2`)
dev.off()

png("recluster_statsplot_OG_3.png", width = 470, height = 800)
plot(statplot$`3`)
dev.off()

png("recluster_statsplot_OG_4.png", width = 470, height = 800)
plot(statplot$`4`)
dev.off()

png("recluster_statsplot_OG_5.png", width = 470, height = 800)
plot(statplot$`5`)
dev.off()

png("recluster_statsplot_OG_6.png", width = 470, height = 800)
plot(statplot$`6`)
dev.off()

png("recluster_statsplot_OG_7.png", width = 470, height = 800)
plot(statplot$`7`)
dev.off()

png("recluster_statsplot_OG_8.png", width = 470, height = 800)
plot(statplot$`8`)
dev.off()

png("recluster_statsplot_OG_9.png", width = 470, height = 800)
plot(statplot$`9`)
dev.off()

png("recluster_statsplot_OG_10.png", width = 470, height = 800)
plot(statplot$`10`)
dev.off()

png("recluster_statsplot_OG_11.png", width = 470, height = 800)
plot(statplot$`11`)
dev.off()

png("recluster_statsplot_OG_12.png", width = 470, height = 800)
plot(statplot$`12`)
dev.off()

png("recluster_statsplot_OG_13.png", width = 470, height = 800)
plot(statplot$`13`)
dev.off()

png("recluster_statsplot_OG_14.png", width = 470, height = 800)
plot(statplot$`14`)
dev.off()

png("recluster_statsplot_OG_15.png", width = 470, height = 800)
plot(statplot$`15`)
dev.off()

png("recluster_statsplot_OG_16.png", width = 470, height = 800)
plot(statplot$`16`)
dev.off()

png("recluster_statsplot_OG_17.png", width = 470, height = 800)
plot(statplot$`17`)
dev.off()

png("recluster_statsplot_OG_18.png", width = 470, height = 800)
plot(statplot$`18`)
dev.off()

## stats that write to table ##
num <- 1:71

## loop to write stats to table 
names <- list()
ids <- list()
zero <- list()
one <- list()
two <- list()
FvalANOVA <- list()
PvalANOVA <- list()
Tukey_Pval_zero_one <- list()
Tukey_Pval_one_two <- list()
Tukey_Pval_two_zero <- list()
KruskalWallis <- list()
PvalKruskal <- list()
Wilcox_Pval_zero_one <- list()
Wilcox_Pval_two_zero <- list()
Wilcox_Pval_two_one <- list()

for (i in num){
  
  n <- eso_proteins$name[i]
  names <- append(names, n)
  id <- eso_proteins$ID[i]
  ids <- append(ids, id)
  
  newdat2 <- dat2[id,]
  
  cn2 <- colnames(newdat2)
  
  setDT(newdat2, keep.rownames = TRUE)[]
  
  ogcountdf <- pivot_longer(newdat2, cols = cn2, names_to = "cell")
  
  colnames(ogcountdf) <- c('ID', 'cell', 'expression')
  
  ogdata <- merge(ogcountdf, ogdf, by = 'cell')
  
  #Mean expression for each group:
  m <- ogdata %>% 
        group_by(cluster) %>% 
        summarise(mean(expression))
  
  z <- m$`mean(expression)`[1]
  zero <- append(zero, z)
  o <- m$`mean(expression)`[2]
  one <- append(one, o)
  t <- m$`mean(expression)`[3]
  two <- append(two, t)
  
  #Anova to compare expression across the three clusters:
  mod <- aov(expression ~ cluster, data = ogdata)
  #View anova:
  sum <- summary(mod)
  f <- sum[[1]][["F value"]][1]
  FvalANOVA <- append(f, FvalANOVA)
  p <- sum[[1]][["Pr(>F)"]][1]
  PvalANOVA <- append(p, PvalANOVA)
  
  tukey <- TukeyHSD(mod)
  tukeydf <- data.frame(tukey$cluster)
  z_o <- tukeydf$p.adj[1]
  Tukey_Pval_zero_one<- append(Tukey_Pval_zero_one, z_o)
  z_t <- tukeydf$p.adj[2]
  Tukey_Pval_one_two <- append(Tukey_Pval_zero_two, z_t)
  t_o <- tukeydf$p.adj[3]
  Tukey_Pval_two_zero <- append(Tukey_Pval_two_one, t_o)
  
  k <- kruskal.test(expression ~ cluster, data = ogdata)
  chi2 <- k[["statistic"]][["Kruskal-Wallis chi-squared"]]
  KruskalWallis <- append(KruskalWallis, chi2)
  pkw <- k[["p.value"]]
  PvalKruskal <- append(PvalKruskal, pkw)
  
  w <- pairwise.wilcox.test(ogdata$expression, ogdata$cluster,
                             p.adjust.method = "none", paired = FALSE)
  wdf <- data.frame(w$p.value)
  w_zo <- wdf$X0[1]
  Wilcox_Pval_zero_one <- append(Wilcox_Pval_zero_one, w_zo)
  w_ot <- wdf$X1[2]
  Wilcox_Pval_one_two <- append(Wilcox_Pval_one_two, w_ot)
  w_tz <- wdf$X0[2]
  Wilcox_Pval_two_zero <- append(Wilcox_Pval_two_zero, w_tz)

}

zero <- as.numeric(zero)
zero <- round(zero, digits = 2)

one <- as.numeric(one)
one <- round(one, digits = 2)

two <- as.numeric(two)
two <- round(two, digits = 2)

FvalANOVA<- as.numeric(FvalANOVA)
FvalANOVA <- round(FvalANOVA, digits = 0)

PvalANOVA<- as.numeric(PvalANOVA)
PvalANOVA<- round(PvalANOVA, digits = 4)

Tukey_Pval_zero_one<- as.numeric(Tukey_Pval_zero_one)
Tukey_Pval_zero_one<- round(Tukey_Pval_zero_one, digits = 4)

Tukey_Pval_two_one<- as.numeric(Tukey_Pval_two_one)
Tukey_Pval_two_one<- round(Tukey_Pval_one_two, digits = 4)

Tukey_Pval_two_zero<- as.numeric(Tukey_Pval_two_zero)
Tukey_Pval_two_zero<- round(Tukey_Pval_two_zero, digits = 4)

KruskalWallis<- as.numeric(KruskalWallis)
KruskalWallis<- round(KruskalWallis, digits = 0)

PvalKruskal<- as.numeric(PvalKruskal)
PvalKruskal<- round(PvalKruskal, digits = 4)

Wilcox_Pval_zero_one<- as.numeric(Wilcox_Pval_zero_one)
Wilcox_Pval_zero_one<- round(Wilcox_Pval_zero_one, digits = 4)

Wilcox_Pval_one_two<- as.numeric(Wilcox_Pval_one_two)
Wilcox_Pval_one_two<- round(Wilcox_Pval_one_two, digits = 4)

Wilcox_Pval_two_zero<- as.numeric(Wilcox_Pval_two_zero)
Wilcox_Pval_two_zero<- round(Wilcox_Pval_two_zero, digits = 4)

s_df3 <- rbind(names,
              ids,
              zero,
              one,
              two,
              FvalANOVA,
              PvalANOVA,
              Tukey_Pval_zero_one,
              Tukey_Pval_two_zero,
              Tukey_Pval_two_one,
              KruskalWallis, 
              PvalKruskal, 
              Wilcox_Pval_zero_one, 
              Wilcox_Pval_one_two, 
              Wilcox_Pval_two_zero)

s_df3 <- data.frame(s_df3)
s_df3<- rotate_df(s_df3)

#write_xlsx(s_df3, "OGrecluster_stats_table.xlsx")



######
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

#3. Plot markers on whole cell umap
num <- 1:26
all.genes <- rownames(seurat_object)
seurat_sub <- ScaleData(seurat_object, features = all.genes)

for (i in num){
  png(sprintf("%s.png", paste('markervln', i)), width = 850, height = 500)
  plot(VlnPlot(seurat_object, features = markersannotatedOG$gene[i], pt.size = 0) +
                                            theme(legend.position = "none", axis.text.x = element_text(angle = 90), plot.title = element_text(size=10))+
                                            labs(title = paste(markersannotatedOG$name[i],"\n", markersannotatedOG$gene[i])))
  dev.off()
}


#need to change for slot thats not scale data
# for (i in num){
#   
#   plot(FeaturePlot(seurat_object, features = markersannotatedOG$gene[i], pt.size = 0.5) +
#                                               scale_colour_gradientn(colours = c("#ABD9E9", "#E0F3F8", "#FDAE61", "#F46D43", "#D73027"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), name = paste('Expression')) +
#                                               theme(plot.title = element_text(size=10))+
#                                               labs(title = paste(markersannotatedOG$name[i],"\n", markersannotatedOG$gene[i])))
#   
#   dev.off()
# }


#SAVE PLOTS

genes <- c('Smp-051920','Smp-105360','Smp-175590','Smp-125320','Smp-010550', 'Smp-172180', 'Smp-085840', 'Smp-124000')
names <- c('nanos-2','notch','fgfra','meg-9','meg-15', 'meg-8.2', 'meg-4.2','meg-14')
stat <- cbind(genes, names)
stat <- data.frame(stat)
num <- 1:8
#550 550
for (i in num){
  png(sprintf("%s.png", paste('reclusterOG', i)), width = 450, height = 450)
  plot(FeaturePlot(seurat_sub, features = stat$genes[i],  slot = 'data', pt.size = 2) + 
    labs(title = paste(stat$names[i],"\n", stat$genes[i]) )+
    theme(plot.title = element_text(size=10)))
  dev.off()
} 

for (i in num){
  png(sprintf("%s.png", paste('reclusterOGvln', i)), width = 850, height = 500)
  plot(VlnPlot(seurat_sub, features = stat$genes[i], pt.size = 0.2, cols = c('#7AD169', '#928DE0', '#FAA637')) +
         theme(legend.position = "none", plot.title = element_text(size=10))+
         labs(title = paste(stat$names[i],"\n", stat$genes[i])))
  dev.off()
}

markersannotatedOG$name
num <- 1:26
for (i in num){
  png(sprintf("%s.png", paste('OGmarker', i)), width = 530, height = 430)
  plot(FeaturePlot(seurat_object, features = markersannotatedOG$gene[i], slot = 'scale.data') + 
         scale_colour_gradientn(colours = c("#74acd1", "#dff3f8", "#fddf90", "#fcad60", "#F46D43", "#d62f27", "#a40025"), breaks=c(-2,0,2,4,6,8,10),labels=c(-2,0,2,4,6,8,10),limits=c(-2,10), na.value = 'white', name = paste('Expression', '\n', 'Z-score')) +
         labs(title = paste(markersannotatedOG$name[i],"\n", markersannotatedOG$gene[i]))+
         theme(plot.title = element_text(size=11)))
  dev.off()
}




