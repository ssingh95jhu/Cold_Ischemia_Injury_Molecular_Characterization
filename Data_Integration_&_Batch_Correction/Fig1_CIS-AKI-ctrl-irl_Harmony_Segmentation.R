#Clean environment
rm(list=ls(all.names=TRUE))
gc()

#Set the working directory
setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

## load data
load("data/AKI_data.RData")
load("data/CIS_data.RData")
load("data/Rabb_ctrl.RData")
load("data/Rabb_irl24h.RData")

#Loading the required library
library(STdeconvolve) #For finding out overdispersed genes
library(harmony) #For implementing harmony to correct for batch effects
library(MUDAN) #To normalize variance?
#library(monocle3) #For trajectory analysis


################################################################################
################# COMBINED CIS,AKI, CTRL AND IRL DATASET SEGMENTATION ##########

##Removing all mitochondrial genes from the gene experssion table (AKI dataset)
AKI_sham$gexp<-AKI_sham$gexp[-grep("^mt", rownames(AKI_sham$gexp)),]
AKI_4h$gexp<-AKI_4h$gexp[-grep("^mt", rownames(AKI_4h$gexp)),]
AKI_12h$gexp<-AKI_12h$gexp[-grep("^mt", rownames(AKI_12h$gexp)),]
AKI_2d$gexp<-AKI_2d$gexp[-grep("^mt", rownames(AKI_2d$gexp)),]
AKI_6w$gexp<-AKI_6w$gexp[-grep("^mt", rownames(AKI_6w$gexp)),]

##Removing all mitochondrial genes from the gene experssion table (irl dataset)
irl1$gexp<-irl1$gexp[-grep("^mt", rownames(irl1$gexp)),]
irl2$gexp<-irl2$gexp[-grep("^mt", rownames(irl2$gexp)),]
irl3$gexp<-irl3$gexp[-grep("^mt", rownames(irl3$gexp)),]
irl4$gexp<-irl4$gexp[-grep("^mt", rownames(irl4$gexp)),]

##Removing all mitochondrial genes from the gene experssion table (ctrl dataset)
ctrl1$gexp<-ctrl1$gexp[-grep("^mt", rownames(ctrl1$gexp)),]
ctrl2$gexp<-ctrl2$gexp[-grep("^mt", rownames(ctrl2$gexp)),]
ctrl3$gexp<-ctrl3$gexp[-grep("^mt", rownames(ctrl3$gexp)),]
ctrl4$gexp<-ctrl4$gexp[-grep("^mt", rownames(ctrl4$gexp)),]

#CIS dataset doesn't have any mitochondrial genes

################################################################################

#Finding overdispersed genes in individual samples

#Overdispersed genes in CIS 0h group
counts <- cleanCounts(CIS_0h$gexp, min.lib.size = 100, min.reads = 10)
corpus_CIS0h <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA, plot=TRUE)
dim(corpus_CIS0h)

#Overdispersed genes in CIS 12h group
counts <- cleanCounts(CIS_12h$gexp, min.lib.size = 100, min.reads = 10)
corpus_CIS12h <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_CIS12h)

#Overdispersed genes in CIS 24h group
counts <- cleanCounts(CIS_24h$gexp, min.lib.size = 100, min.reads = 10)
corpus_CIS24h <- restrictCorpus(counts, removeAbove=1.1, removeBelow = 0.05, nTopOD = NA)
dim(corpus_CIS24h)

#Overdispersed genes in CIS 48h group
counts <- cleanCounts(CIS_48h$gexp, min.lib.size = 100, min.reads = 10)
corpus_CIS48h <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_CIS48h)

#Overdispersed genes in AKI sham group
counts <- cleanCounts(AKI_sham$gexp, min.lib.size = 100, min.reads = 10)
corpus_AKIsham <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_AKIsham)

#Overdispersed genes in AKI 4h group
counts <- cleanCounts(AKI_4h$gexp, min.lib.size = 100, min.reads = 10)
corpus_AKI4h <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_AKI4h)

#Overdispersed genes in AKI 12h group
counts <- cleanCounts(AKI_12h$gexp, min.lib.size = 100, min.reads = 10)
corpus_AKI12h <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_AKI12h)

#Overdispersed genes in AKI 2d group
counts <- cleanCounts(AKI_2d$gexp, min.lib.size = 100, min.reads = 10)
corpus_AKI2d <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_AKI2d)

#Overdispersed genes in AKI 6w group
counts <- cleanCounts(AKI_6w$gexp, min.lib.size = 100, min.reads = 10)
corpus_AKI6w <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_AKI6w)

#Overdispersed genes in irl1 group
counts <- cleanCounts(irl1$gexp, min.lib.size = 100, min.reads = 10)
corpus_irl1 <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_irl1)

#Overdispersed genes in irl2 group
counts <- cleanCounts(irl2$gexp, min.lib.size = 100, min.reads = 10)
corpus_irl2 <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_irl2)


#Overdispersed genes in irl3 group
counts <- cleanCounts(irl3$gexp, min.lib.size = 100, min.reads = 10)
corpus_irl3 <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_irl3)

#Overdispersed genes in irl4 group
counts <- cleanCounts(irl4$gexp, min.lib.size = 100, min.reads = 10)
corpus_irl4 <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_irl4)

#Overdispersed genes in ctrl1 group
counts <- cleanCounts(ctrl1$gexp, min.lib.size = 100, min.reads = 10)
corpus_ctrl1 <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_ctrl1)

#Overdispersed genes in ctrl2 group
counts <- cleanCounts(ctrl2$gexp, min.lib.size = 100, min.reads = 10)
corpus_ctrl2 <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_ctrl2)

#Overdispersed genes in ctrl3 group
counts <- cleanCounts(ctrl3$gexp, min.lib.size = 100, min.reads = 10)
corpus_ctrl3 <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_ctrl3)

#Overdispersed genes in ctrl4 group.
counts <- cleanCounts(ctrl4$gexp, min.lib.size = 100, min.reads = 10)
corpus_ctrl4 <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = NA)
dim(corpus_ctrl4)

data<-c(rownames(corpus_CIS0h),rownames(corpus_CIS12h), 
        rownames(corpus_CIS24h), rownames(corpus_CIS48h),
        rownames(corpus_AKIsham),rownames(corpus_AKI4h),rownames(corpus_AKI12h), 
        rownames(corpus_AKI2d), rownames(corpus_AKI6w),
        rownames(corpus_irl1),rownames(corpus_irl2), 
        rownames(corpus_irl3), rownames(corpus_irl4),
        rownames(corpus_ctrl1),rownames(corpus_ctrl2), 
        rownames(corpus_ctrl3),rownames(corpus_ctrl4) )

#COD<-data[which(table(data)>=2)] <- not the desired set of gene names
COD<-unique(names(which(table(data)>=2)))
COD<-intersect(COD,rownames(CIS_0h$gexp))
#saveRDS(COD, file= "OverDispersed_Feautures_for_PC.rds")

library(openxlsx)
wb1<-createWorkbook()
addWorksheet(wb1,"Overdispersed_genes")
writeData(wb1, "Overdispersed_genes",COD)
saveWorkbook(wb1,"CIS_AKI_CTRL_IRL_Top_Feature_Selected_Genes.xlsx")

## Impelementing HARMONY
#Getting the data
dataA<-CIS_0h$gexp[COD,]
dataB<-CIS_12h$gexp[COD,]
dataC<-CIS_24h$gexp[COD,]
dataD<-CIS_48h$gexp[COD,]

dataE<-AKI_sham$gexp[COD,]
dataF<-AKI_4h$gexp[COD,]
dataG<-AKI_12h$gexp[COD,]
dataH<-AKI_2d$gexp[COD,]
dataI<-AKI_6w$gexp[COD,]

dataJ<-irl1$gexp[COD,]
dataK<-irl2$gexp[COD,]
dataL<-irl3$gexp[COD,]
dataM<-irl4$gexp[COD,]

dataN<-ctrl1$gexp[COD,]
dataO<-ctrl2$gexp[COD,]
dataP<-ctrl3$gexp[COD,]
dataQ<-ctrl4$gexp[COD,]

#Combine different datasets into a single dataset
cd <- cbind(dataA, dataB, dataC, dataD, 
            dataE, dataF, dataG, dataH, dataI,
            dataJ, dataK, dataL, dataM,
            dataN, dataO, dataP, dataQ)
cd<-as.matrix(cd)
#cd<-cd[,-c(19931,20077)]
cd<-cd[,-c(19931,20077)] #These two spots have no gene expression in them.

# meta data 
meta <- c( rep('CIS 0h', ncol(dataA)), rep('CIS 12h', ncol(dataB)), 
           rep('CIS_24h', ncol(dataC)), rep('CIS_48h', ncol(dataD)), 
           rep('AKI_sham', ncol(dataE)), rep('AKI_4h', ncol(dataF)), rep('AKI_12h', ncol(dataG)),
           rep('AKI_2d', ncol(dataH)), rep('AKI_6w', ncol(dataI)),
          
           rep('irl1', ncol(dataJ)), rep('irl2', ncol(dataK)), 
           rep('irl3', ncol(dataL)), rep('irl4', ncol(dataM)),
           rep('ctrl1', ncol(dataN)), rep('ctrl2', ncol(dataO)), 
           rep('ctrl3', ncol(dataP)), rep('ctrl4', ncol(dataQ)))

names(meta) <- c(colnames(dataA), colnames(dataB), colnames(dataC), colnames(dataD),
                 colnames(dataE), colnames(dataF), colnames(dataG), colnames(dataH), colnames(dataI),
                 colnames(dataJ), colnames(dataK), colnames(dataL), colnames(dataM),
                 colnames(dataN), colnames(dataO), colnames(dataP), colnames(dataQ))
meta<-meta[-c(19931,20077)]

meta <- factor(meta)


#mat <- MERINGUE::cleanCounts(cd)
#dim(mat)
#dim(cd)
## CPM normalization (is it the same as in MUDAN?)
matnorm <- MERINGUE::normalizeCounts(cd, verbose=FALSE) 
table(is.na(as.vector(matnorm)))

pcs <- MUDAN::getPcs(matnorm,
                     nGenes=length(matnorm), 
                     nPcs=30,
                     verbose=FALSE) 
#write.csv(pcs,file= "AKI-CIS-irl-ctrl_pcs.csv")
saveRDS(pcs, file = "AKI-CIS-irl-ctrl_pcs.rds")

#pcs<-prcomp(na.omit(matnorm), scale=F)


#UMAP Plotting of regular (non-harmonized) clusters
#emb <- uwot::umap(pcs) ## 2d embedding
#rownames(emb) <- colnames(matnorm)

emb<- Rtsne::Rtsne(pcs)
emb <- emb$Y
rownames(emb) <- colnames(matnorm)
head(emb)
saveRDS(emb, file = "AKI-CIS-irl-ctrl_emb.rds")

pdf('Non-Harmonized_AKI-CIS-IRL-CTRL_datasets.pdf', height=5, width=5)
par(mfrow=c(1,1))
plotEmbedding(emb, groups=meta, 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              #main='Non-harmonized CIS+irl+ctrl dataset',
              verbose=FALSE,legend.x='topleft', 
              mark.clusters = TRUE, mark.cluster.cex = 0.6,
              shuffle.colors = FALSE, cex=0.3, alpha=1)
dev.off()


# Now harmonize PCs
harmonized <- HarmonyMatrix(pcs, meta, do_pca = FALSE, verbose = FALSE,
                            theta = 10, lambda = 0.05)
saveRDS(harmonized, file = "AKI-CIS-irl-ctrl_Harmonized_Matrix.rds")

emb.harmony <- Rtsne::Rtsne(harmonized)
emb.harmony <- emb.harmony$Y
rownames(emb.harmony) <- colnames(matnorm)
head(emb.harmony) 


#write.csv(emb.harmony,file= "AKI-CIS-irl-ctrl_Harmonized_tSNE_Embeds.csv")
saveRDS(emb.harmony, file = "AKI-CIS-irl-ctrl_Harmonized_tSNE_Embeds.rds")

pdf('Harmonized_AKI-CIS-IRL-CTRL_datasets.pdf', height=5, width=5)
plotEmbedding(emb.harmony, groups=meta, 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              #main='Non-harmonized CIS+irl+ctrl dataset',
              verbose=FALSE,legend.x='topleft', 
              mark.clusters = FALSE, mark.cluster.cex = 1,
              shuffle.colors = FALSE, cex=0.3, alpha=1)
dev.off()

# Joint clustering
#com <- getComMembership(harmonized, k=100, method=igraph::cluster_walktrap)
#com <- getComMembership(harmonized, k=100, method=igraph::louvain)
com <- getApproxComMembership(harmonized, k=100, method=igraph::cluster_louvain)

#Saving the file
#write.csv(com, file= "AKI-CIS-irl-ctrl_Harmonized_AllClusters.csv")
saveRDS(com, file = "AKI-CIS-irl-ctrl_Harmonized_AllClusters.rds")

##################################################################
# com<-read.csv(file= "AKI-CIS-irl-ctrl_Harmonized_AllClusters.csv")
# rownames(com)<-com$X
# com$X<-com$x
# head(com)
# com<-subset(com, select = -x)
# head(com)
# com<-as.matrix(com)
# com<-as.factor(com)
###################################################################

pdf('Harmonized AKI-CIS-IRL-CTRL_TransClusters.pdf', height=5, width=5)
par(mfrow=c(1,1), mar=rep(2,4))
plotEmbedding(emb.harmony, groups=com, 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              #main='Harmonized CIS+irl+ctrl dataset',
              verbose=FALSE,legend.x='topleft', 
              mark.clusters = TRUE, mark.cluster.cex = 1,
              shuffle.colors = FALSE, cex=0.3, alpha=1)
dev.off()

plotEmbedding(emb.harmony, groups=meta, 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='Harmonized AKI-CIS-irl-ctrl dataset',
              verbose=FALSE,legend.x='topleft', 
              #mark.clusters = TRUE, mark.cluster.cex = 1,
              #shuffle.colors = TRUE, cex=0.1
              )

n<-0.7
#Visualizing the clusters in a spatial context
pdf('CIS_Whole_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(CIS_0h$pos, groups=com[rownames(CIS_0h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 0h', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_12h$pos, groups=com[rownames(CIS_12h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 12h', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_24h$pos, groups=com[rownames(CIS_24h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 24h', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_48h$pos, groups=com[rownames(CIS_48h$pos)], 
              show.legend=TRUE, xlab=NA, ylab=NA, 
              main='CIS 48h', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the clusters in a spatial context (irl)
pdf('IRL_Whole_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(irl1$pos, groups=com[rownames(irl1$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl1', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl2$pos, groups=com[rownames(irl2$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl2', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl3$pos, groups=com[rownames(irl3$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl3', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl4$pos, groups=com[rownames(irl4$pos)], 
              show.legend=TRUE, xlab=NA, ylab=NA, 
              main='irl4', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the clusters in a spatial context (ctrl)
pdf('CTRL_Whole_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(ctrl1$pos, groups=com[rownames(ctrl1$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl1', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl2$pos, groups=com[rownames(ctrl2$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl2', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl3$pos, groups=com[rownames(ctrl3$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl3', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl4$pos, groups=com[rownames(ctrl4$pos)], 
              show.legend=TRUE, xlab=NA, ylab=NA, 
              main='ctrl4', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the clusters in a spatial context (AKI)
pdf('AKI_Whole_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,5), mar=rep(2,4))
plotEmbedding(AKI_sham$pos, groups=com[rownames(AKI_sham$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI sham', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_4h$pos, groups=com[rownames(AKI_4h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 4h', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_12h$pos, groups=com[rownames(AKI_12h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 12h', cex=n, alpha=1,
             verbose=FALSE)

plotEmbedding(AKI_2d$pos, groups=com[rownames(AKI_2d$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 2days', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_6w$pos, groups=com[rownames(AKI_6w$pos)], 
              show.legend=TRUE, xlab=NA, ylab=NA, 
              main='AKI 6 weeks', cex=n, alpha=1,
              verbose=FALSE)
dev.off()


#Segmenting out the medulla region from all the timepoints
medulla<-com[com==1]
saveRDS(medulla, file = "AKI-CIS-irl-ctrl_Medulla_spots.rds")

#Visualizing the segmented medulla. (CIS)
pdf('CIS_Medulla_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(CIS_0h$pos, groups=medulla[rownames(CIS_0h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 0h medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_12h$pos, groups=medulla[rownames(CIS_12h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 12h medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_24h$pos, groups=medulla[rownames(CIS_24h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 24h medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_48h$pos, groups=medulla[rownames(CIS_48h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 48h medulla', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the segmented medulla (irl)
pdf('IRL_Medulla_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(irl1$pos, groups=medulla[rownames(irl1$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl1 medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl2$pos, groups=medulla[rownames(irl2$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl2 medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl3$pos, groups=medulla[rownames(irl3$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl3 medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl4$pos, groups=medulla[rownames(irl4$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl4 medulla', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the segmented medulla (ctrl)
pdf('CTRL_Medulla_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(ctrl1$pos, groups=medulla[rownames(ctrl1$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl1 medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl2$pos, groups=medulla[rownames(ctrl2$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl2 medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl3$pos, groups=medulla[rownames(ctrl3$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl3 medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl4$pos, groups=medulla[rownames(ctrl4$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl4 medulla', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the segmented medulla (aki)
pdf('AKI_Medulla_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,5), mar=rep(2,4))
plotEmbedding(AKI_sham$pos, groups=medulla[rownames(AKI_sham$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI sham medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_4h$pos, groups=medulla[rownames(AKI_4h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 4h medulla ', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_12h$pos, groups=medulla[rownames(AKI_12h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 12h medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_2d$pos, groups=medulla[rownames(AKI_2d$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 2d medulla', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_6w$pos, groups=medulla[rownames(AKI_6w$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 6w medulla', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Segmenting out the interface
interface<-com[com==7|com==8]
saveRDS(interface, file = "AKI-CIS-irl-ctrl_Interface_spots.rds")

#Visualizing the segmented interface (CIS)
pdf('CIS_Interface_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(CIS_0h$pos, groups=interface[rownames(CIS_0h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 0h interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_12h$pos, groups=interface[rownames(CIS_12h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 12h interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_24h$pos, groups=interface[rownames(CIS_24h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 24h interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_48h$pos, groups=interface[rownames(CIS_48h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 48h interface', cex=n, alpha=1,
              verbose=FALSE)
dev.off()


#Visualizing the segmented interface (irl)
pdf('IRL_Interface_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(irl1$pos, groups=interface[rownames(irl1$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl1 interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl2$pos, groups=interface[rownames(irl2$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl2 interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl3$pos, groups=interface[rownames(irl3$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl3 interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl4$pos, groups=interface[rownames(irl4$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl4 interface', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the segmented interface (ctrl)
pdf('CTRL_Interface_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(ctrl1$pos, groups=interface[rownames(ctrl1$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl1 interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl2$pos, groups=interface[rownames(ctrl2$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl2 interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl3$pos, groups=interface[rownames(ctrl3$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl3 interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl4$pos, groups=interface[rownames(ctrl4$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl4 interface', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the segmented interface (AKI)
pdf('AKI_Interface_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,5), mar=rep(2,4))
plotEmbedding(AKI_sham$pos, groups=interface[rownames(AKI_sham$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI sham interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_4h$pos, groups=interface[rownames(AKI_4h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 4h interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_12h$pos, groups=interface[rownames(AKI_12h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 12h interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_2d$pos, groups=interface[rownames(AKI_2d$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 2d interface', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_6w$pos, groups=interface[rownames(AKI_6w$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 6w interface', cex=n, alpha=1,
              verbose=FALSE)
dev.off()


#Segmenting out the cortex
#cortex<-com[com!=10 & com!=6 & com!=1 & com!=8 & com!=9]
cortex<-com[com==2|com==3|com==4|com==6| com==10]
saveRDS(cortex, file = "AKI-CIS-irl-ctrl_Cortex_spots.rds")

#Visualizing the segmented cortex (CIS)
pdf('CIS_Cortex_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(CIS_0h$pos, groups=cortex[rownames(CIS_0h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 0h cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_12h$pos, groups=cortex[rownames(CIS_12h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 12h cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_24h$pos, groups=cortex[rownames(CIS_24h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 24h cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(CIS_48h$pos, groups=cortex[rownames(CIS_48h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='CIS 48h cortex', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the segmented cortex (irl)
pdf('IRL_Cortex_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(irl1$pos, groups=cortex[rownames(irl1$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl1 cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl2$pos, groups=cortex[rownames(irl2$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl2 cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl3$pos, groups=cortex[rownames(irl3$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl3 cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(irl4$pos, groups=cortex[rownames(irl4$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='irl4 cortex', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the segmented cortex (ctrl)
pdf('CTRL_Cortex_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,4), mar=rep(2,4))
plotEmbedding(ctrl1$pos, groups=cortex[rownames(ctrl1$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl1 cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl2$pos, groups=cortex[rownames(ctrl2$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl2 cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl3$pos, groups=cortex[rownames(ctrl3$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl3 cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(ctrl4$pos, groups=cortex[rownames(ctrl4$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='ctrl4 cortex', cex=n, alpha=1,
              verbose=FALSE)
dev.off()

#Visualizing the segmented cortex (aki)
pdf('AKI_Cortex_All_timepoints.pdf',height=3, width=15)
par(mfrow=c(1,5), mar=rep(2,4))
plotEmbedding(AKI_sham$pos, groups=cortex[rownames(AKI_sham$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI sham cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_4h$pos, groups=cortex[rownames(AKI_4h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 4h cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_12h$pos, groups=cortex[rownames(AKI_12h$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 12h cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_2d$pos, groups=cortex[rownames(AKI_2d$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 2d cortex', cex=n, alpha=1,
              verbose=FALSE)

plotEmbedding(AKI_6w$pos, groups=cortex[rownames(AKI_6w$pos)], 
              show.legend=FALSE, xlab=NA, ylab=NA, 
              main='AKI 6w cortex', cex=n, alpha=1,
              verbose=FALSE)
dev.off()
