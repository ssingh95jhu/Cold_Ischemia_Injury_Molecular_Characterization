#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

## load data
load("data/CIS_data.RData")
load("data/AKI_data.RData")
load("data/Rabb_ctrl.RData")

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(msigdbr) #For HALLMARK Pathway 
library(enrichplot) #For ClusterProfiler()
library(pathview)
library(biomaRt)

##############################################################################
#Extracting the compartment specific spots for CIS and AKI:

cortex <- readRDS('AKI-CIS-irl-ctrl_Cortex_spots.rds')
interface <- readRDS('AKI-CIS-irl-ctrl_Interface_spots.rds')
medulla <- readRDS('AKI-CIS-irl-ctrl_Medulla_spots.rds')
#all_spots<-readRDS('AKI-CIS-irl-ctrl_All_spots.rds')

############################## CIS DATASET SPOTS ###############################
#1.All_Cortex_Spots
CIS_cortex.0h<-cortex[grep("^CIS_0h", names(cortex))]
CIS_cortex.12h<-cortex[grep("^CIS_12h", names(cortex))]
CIS_cortex.24h<-cortex[grep("^CIS_24h", names(cortex))]
CIS_cortex.48h<-cortex[grep("^CIS_48h", names(cortex))]

CIS_cortex.0h<-cbind(CIS_0h$pos[names(CIS_cortex.0h),], CIS_cortex.0h)
CIS_cortex.12h<-cbind(CIS_12h$pos[names(CIS_cortex.12h),], CIS_cortex.12h)
CIS_cortex.24h<-cbind(CIS_24h$pos[names(CIS_cortex.24h),], CIS_cortex.24h)
CIS_cortex.48h<-cbind(CIS_48h$pos[names(CIS_cortex.48h),], CIS_cortex.48h)

colnames(CIS_cortex.0h)<-c('X','Y', 'Cluster')
colnames(CIS_cortex.12h)<-c('X','Y', 'Cluster')
colnames(CIS_cortex.24h)<-c('X','Y', 'Cluster')
colnames(CIS_cortex.48h)<-c('X','Y', 'Cluster')

#To plot and check/verify the spots are right
ggplot(data.frame(CIS_cortex.0h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_cortex_0h")
ggplot(data.frame(CIS_cortex.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_cortex_12h")
ggplot(data.frame(CIS_cortex.24h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_cortex_24h")
ggplot(data.frame(CIS_cortex.48h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_cortex_48h")

#2.All Interface spots
CIS_interface.0h<-interface[grep("^CIS_0h", names(interface))]
CIS_interface.12h<-interface[grep("^CIS_12h", names(interface))]
CIS_interface.24h<-interface[grep("^CIS_24h", names(interface))]
CIS_interface.48h<-interface[grep("^CIS_48h", names(interface))]

CIS_interface.0h<-cbind(CIS_0h$pos[names(CIS_interface.0h),], CIS_interface.0h)
CIS_interface.12h<-cbind(CIS_12h$pos[names(CIS_interface.12h),], CIS_interface.12h)
CIS_interface.24h<-cbind(CIS_24h$pos[names(CIS_interface.24h),], CIS_interface.24h)
CIS_interface.48h<-cbind(CIS_48h$pos[names(CIS_interface.48h),], CIS_interface.48h)

colnames(CIS_interface.0h)<-c('X','Y', 'Cluster')
colnames(CIS_interface.12h)<-c('X','Y', 'Cluster')
colnames(CIS_interface.24h)<-c('X','Y', 'Cluster')
colnames(CIS_interface.48h)<-c('X','Y', 'Cluster')

#To plot and check/verify the spots are right
ggplot(data.frame(CIS_interface.0h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_Interface_0h")
ggplot(data.frame(CIS_interface.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_Interface_12h")
ggplot(data.frame(CIS_interface.24h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_Interface_24h")
ggplot(data.frame(CIS_interface.48h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_Interface_48h")

#3.All medullary spots
CIS_medulla.0h<-medulla[grep("^CIS_0h", names(medulla))]
CIS_medulla.12h<-medulla[grep("^CIS_12h", names(medulla))]
CIS_medulla.24h<-medulla[grep("^CIS_24h", names(medulla))]
CIS_medulla.48h<-medulla[grep("^CIS_48h", names(medulla))]

CIS_medulla.0h<-cbind(CIS_0h$pos[names(CIS_medulla.0h),], CIS_medulla.0h)
CIS_medulla.12h<-cbind(CIS_12h$pos[names(CIS_medulla.12h),], CIS_medulla.12h)
CIS_medulla.24h<-cbind(CIS_24h$pos[names(CIS_medulla.24h),], CIS_medulla.24h)
CIS_medulla.48h<-cbind(CIS_48h$pos[names(CIS_medulla.48h),], CIS_medulla.48h)

colnames(CIS_medulla.0h)<-c('X','Y', 'Cluster')
colnames(CIS_medulla.12h)<-c('X','Y', 'Cluster')
colnames(CIS_medulla.24h)<-c('X','Y', 'Cluster')
colnames(CIS_medulla.48h)<-c('X','Y', 'Cluster')

#To plot and check/verify the spots are right
ggplot(data.frame(CIS_medulla.0h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_Medulla_0h")
ggplot(data.frame(CIS_medulla.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_Medulla_12h")
ggplot(data.frame(CIS_medulla.24h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_Medulla_24h")
ggplot(data.frame(CIS_medulla.48h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("CIS_Medulla_48h")

################################################################################
############################ AKI SPOTS #########################################
#1.All_Cortex_Spots
AKI_cortex.sham<-cortex[grep("^AKI_sham", names(cortex))]
AKI_cortex.4h<-cortex[grep("^AKI_4h", names(cortex))]
AKI_cortex.12h<-cortex[grep("^AKI_12h", names(cortex))]
AKI_cortex.2d<-cortex[grep("^AKI_2d", names(cortex))]

AKI_cortex.sham<-cbind(AKI_sham$pos[names(AKI_cortex.sham),], AKI_cortex.sham)
AKI_cortex.4h<-cbind(AKI_4h$pos[names(AKI_cortex.4h),], AKI_cortex.4h)
AKI_cortex.12h<-cbind(AKI_12h$pos[names(AKI_cortex.12h),], AKI_cortex.12h)
AKI_cortex.2d<-cbind(AKI_2d$pos[names(AKI_cortex.2d),], AKI_cortex.2d)

colnames(AKI_cortex.sham)<-c('X','Y', 'Cluster')
colnames(AKI_cortex.4h)<-c('X','Y', 'Cluster')
colnames(AKI_cortex.12h)<-c('X','Y', 'Cluster')
colnames(AKI_cortex.2d)<-c('X','Y', 'Cluster')

#To plot and check/verify the spots are right
ggplot(data.frame(AKI_cortex.sham), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Cortex_sham")
ggplot(data.frame(AKI_cortex.4h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Cortex_4h")
ggplot(data.frame(AKI_cortex.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Cortex_12h")
ggplot(data.frame(AKI_cortex.2d), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Cortex_2d")

#2.All Interface spots
AKI_interface.sham<-interface[grep("^AKI_sham", names(interface))]
AKI_interface.4h<-interface[grep("^AKI_4h", names(interface))]
AKI_interface.12h<-interface[grep("^AKI_12h", names(interface))]
AKI_interface.2d<-interface[grep("^AKI_2d", names(interface))]

AKI_interface.sham<-cbind(AKI_sham$pos[names(AKI_interface.sham),], AKI_interface.sham)
AKI_interface.4h<-cbind(AKI_4h$pos[names(AKI_interface.4h),], AKI_interface.4h)
AKI_interface.12h<-cbind(AKI_12h$pos[names(AKI_interface.12h),], AKI_interface.12h)
AKI_interface.2d<-cbind(AKI_2d$pos[names(AKI_interface.2d),], AKI_interface.2d)

colnames(AKI_interface.sham)<-c('X','Y', 'Cluster')
colnames(AKI_interface.4h)<-c('X','Y', 'Cluster')
colnames(AKI_interface.12h)<-c('X','Y', 'Cluster')
colnames(AKI_interface.2d)<-c('X','Y', 'Cluster')

#To plot and check/verify the spots are right
ggplot(data.frame(AKI_interface.sham), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Interface_sham")
ggplot(data.frame(AKI_interface.4h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Interface_4h")
ggplot(data.frame(AKI_interface.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Interface_12h")
ggplot(data.frame(AKI_interface.2d), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Interface_2d")


#3.All medullary spots
AKI_medulla.sham<-medulla[grep("^AKI_sham", names(medulla))]
AKI_medulla.4h<-medulla[grep("^AKI_4h", names(medulla))]
AKI_medulla.12h<-medulla[grep("^AKI_12h", names(medulla))]
AKI_medulla.2d<-medulla[grep("^AKI_2d", names(medulla))]

AKI_medulla.sham<-cbind(AKI_sham$pos[names(AKI_medulla.sham),], AKI_medulla.sham)
AKI_medulla.4h<-cbind(AKI_4h$pos[names(AKI_medulla.4h),], AKI_medulla.4h)
AKI_medulla.12h<-cbind(AKI_12h$pos[names(AKI_medulla.12h),], AKI_medulla.12h)
AKI_medulla.2d<-cbind(AKI_2d$pos[names(AKI_medulla.2d),], AKI_medulla.2d)

colnames(AKI_medulla.sham)<-c('X','Y', 'Cluster')
colnames(AKI_medulla.4h)<-c('X','Y', 'Cluster')
colnames(AKI_medulla.12h)<-c('X','Y', 'Cluster')
colnames(AKI_medulla.2d)<-c('X','Y', 'Cluster')

#To plot and check/verify the spots are right
ggplot(data.frame(AKI_medulla.sham), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Medulla_sham")
ggplot(data.frame(AKI_medulla.4h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Medulla_4h")
ggplot(data.frame(AKI_medulla.12h), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Medulla_12h")
ggplot(data.frame(AKI_medulla.2d), aes(x=X, y=Y, col=Cluster)) + geom_point(size=0.5) + ggtitle("AKI_Medulla_2d")

################################################################################


#Extracting all CIS spots excluding the ureter
# CIS_allspots.0h<-c(rownames(CIS_cortex.0h),rownames(CIS_interface.0h),rownames(CIS_medulla.0h))
# CIS_allspots.12h<-c(rownames(CIS_cortex.12h),rownames(CIS_interface.12h),rownames(CIS_medulla.12h))
# CIS_allspots.24h<-c(rownames(CIS_cortex.24h),rownames(CIS_interface.24h),rownames(CIS_medulla.24h))
# CIS_allspots.48h<-c(rownames(CIS_cortex.48h),rownames(CIS_interface.48h),rownames(CIS_medulla.48h))


################################################################################
#################### NORMALIZE CTRL DATASET ####################################

## limma performs better with gaussian distributed
ctrl1$mat <- MERINGUE::normalizeCounts(ctrl1$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
ctrl2$mat <- MERINGUE::normalizeCounts(ctrl2$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
ctrl3$mat <- MERINGUE::normalizeCounts(ctrl3$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
ctrl4$mat <- MERINGUE::normalizeCounts(ctrl4$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)

## easier for plotting
ctrl1$mat_notlog <- MERINGUE::normalizeCounts(ctrl1$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl2$mat_notlog <- MERINGUE::normalizeCounts(ctrl2$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl3$mat_notlog <- MERINGUE::normalizeCounts(ctrl3$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl4$mat_notlog <- MERINGUE::normalizeCounts(ctrl4$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

############################ normalize #########################################

## limma performs better with gaussian distributed
CIS_0h$mat <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
CIS_12h$mat <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
CIS_24h$mat <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
CIS_48h$mat <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)

## easier for plotting
CIS_0h$mat_notlog <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_12h$mat_notlog <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_24h$mat_notlog <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_48h$mat_notlog <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

###################### AKI DATASET #############################################

AKI_2d$gexp<-AKI_2d$gexp[,-1818] #This spot has no genes
AKI_2d$pos<-AKI_2d$pos[colnames(AKI_2d$gexp),]

##Removing all mitochondrial genes from the gene experssion table
AKI_sham$gexp<-AKI_sham$gexp[-grep("^mt", rownames(AKI_sham$gexp)),]
AKI_4h$gexp<-AKI_4h$gexp[-grep("^mt", rownames(AKI_4h$gexp)),]
AKI_12h$gexp<-AKI_12h$gexp[-grep("^mt", rownames(AKI_12h$gexp)),]
AKI_2d$gexp<-AKI_2d$gexp[-grep("^mt", rownames(AKI_2d$gexp)),]
AKI_6w$gexp<-AKI_6w$gexp[-grep("^mt", rownames(AKI_6w$gexp)),]

##Updating the position matrix for AKI
AKI_sham$pos<-AKI_sham$pos[grep("^AKI_sham", colnames(AKI_sham$gexp)),]
AKI_4h$pos<-AKI_4h$pos[grep("^AKI_4h", colnames(AKI_4h$gexp)),]
AKI_12h$pos<-AKI_12h$pos[grep("^AKI_12h", colnames(AKI_12h$gexp)),]
AKI_2d$pos<-AKI_2d$pos[grep("^AKI_2d", colnames(AKI_2d$gexp)),]
AKI_6w$pos<-AKI_6w$pos[grep("^AKI_6w", colnames(AKI_6w$gexp)),]

################# normalize AKI DATASET (CPM) ##################################

AKI_sham$mat <- MERINGUE::normalizeCounts(AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
AKI_4h$mat <- MERINGUE::normalizeCounts(AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
AKI_12h$mat <- MERINGUE::normalizeCounts(AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
AKI_2d$mat <- MERINGUE::normalizeCounts(AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
AKI_6w$mat <- MERINGUE::normalizeCounts(AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)

AKI_sham$mat_notlog <- MERINGUE::normalizeCounts(AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_4h$mat_notlog <- MERINGUE::normalizeCounts(AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_12h$mat_notlog <- MERINGUE::normalizeCounts(AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_2d$mat_notlog <- MERINGUE::normalizeCounts(AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_6w$mat_notlog <- MERINGUE::normalizeCounts(AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

################################################################################
# #Creating background gene list for GSEA 
# gene_sets<-msigdbr(species = "mouse", category = "H")
# background <- gene_sets[, c("gs_name", "gene_symbol")]



############### LINEAR REGRESSION ANALYSIS #####################################
######################### GLOBAL DEG (CIS) #####################################
#Global refers to the whole tissue section without any compartmentalization.
#In the main article section, 
#TIME ENCODED AS LOG (HOURS + 1)
library(limma)

wb1<-createWorkbook() #workbook for DEG genelist (CIS 0 - 48hrs)
wb2<-createWorkbook() #workbook for gsea GO pathway (CIS 0 - 48hrs)
wb3<-createWorkbook() #workbook for gsea HALLMARK pathway (CIS 0 - 48 hrs)
wb4<-createWorkbook() #workbook for gsea KEGG pathway (CIS 0 - 48 hrs)

genes.shared <- Reduce(intersect, list(rownames(CIS_0h$mat_notlog), rownames(CIS_12h$mat_notlog),
                                      rownames(CIS_24h$mat_notlog), rownames(CIS_48h$mat_notlog),
                                      rownames(AKI_sham$mat_notlog), rownames(AKI_4h$mat_notlog),
                                      rownames(AKI_12h$mat_notlog), rownames(AKI_2d$mat_notlog)) )
#genes.shared<-rownames(CIS_0h$gexp)

#Selecting n=190 random points from the whole CIS tissue cross section
sample_n<-190 #Setting the no. of random spots to be extracted

set.seed(100)
cis_rand_spots.0h<-sample(colnames(CIS_0h$mat_notlog), sample_n, replace=FALSE)
cis_rand_spots.12h<-sample(colnames(CIS_12h$mat_notlog), sample_n, replace=FALSE)
cis_rand_spots.24h<-sample(colnames(CIS_24h$mat_notlog), sample_n, replace=FALSE)
cis_rand_spots.48h<-sample(colnames(CIS_48h$mat_notlog), sample_n, replace=FALSE)

#Plot to check:
plot(CIS_0h$pos[cis_rand_spots.0h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Global_0h")
plot(CIS_12h$pos[cis_rand_spots.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Global_12h")
plot(CIS_24h$pos[cis_rand_spots.24h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Global_24h")
plot(CIS_48h$pos[cis_rand_spots.48h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Global_48h")


# CIS_gexp <- data.frame(cbind(
#   as.matrix(CIS_0h$mat_notlog[genes.shared,]),
#   as.matrix(CIS_12h$mat_notlog[genes.shared,]),
#   as.matrix(CIS_24h$mat_notlog[genes.shared,]),
#   as.matrix(CIS_48h$mat_notlog[genes.shared,]) ))

CIS_gexp <- data.frame(cbind(
  as.matrix(CIS_0h$mat_notlog[genes.shared,cis_rand_spots.0h]),
  as.matrix(CIS_12h$mat_notlog[genes.shared,cis_rand_spots.12h]),
  as.matrix(CIS_24h$mat_notlog[genes.shared,cis_rand_spots.24h]),
  as.matrix(CIS_48h$mat_notlog[genes.shared,cis_rand_spots.48h]) ))

# meta <- data.frame(time = c(
#   rep(log(0+1), ncol(CIS_0h$mat_notlog)),
#   rep(log(12+1), ncol(CIS_12h$mat_notlog)), ## unit of days
#   rep(log(24+1), ncol(CIS_24h$mat_notlog)),
#   rep(log(48+1), ncol(CIS_48h$mat_notlog))
# ))  
# rownames(meta) <- colnames(CIS_gexp)
# table(meta)

meta <- data.frame(time = c(
  rep(log(0+1), length(cis_rand_spots.0h)),
  rep(log(12+1), length(cis_rand_spots.12h)), ## unit of days
  rep(log(24+1), length(cis_rand_spots.24h)),
  rep(log(48+1), length(cis_rand_spots.24h))
))  
rownames(meta) <- colnames(CIS_gexp)
table(meta)

#des = with(meta, model.matrix(~ time)) 
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(CIS_gexp, design = des) ## fit limma model against time
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(CIS_gexp))
top <- top[top$adj.P.Val < 0.05,] ## restrict to significant results only
head(top)

cis_global_fc <- top$logFC
names(cis_global_fc) <- rownames(top)
cis_global_fc<-data.frame(names(cis_global_fc) ,cis_global_fc)
colnames(cis_global_fc)<-c('Gene','Slope')
head(cis_global_fc)

cis_up<-cis_global_fc[cis_global_fc$Slope>=0,]
head(cis_up)

cis_down<-cis_global_fc[cis_global_fc$Slope<0,]
head(cis_down)

#Creating a ranked list of genes based on their slope
cis_genes<-rbind(cis_up,cis_down)
print(cis_genes)
cis_genes$rank<-c(rank(cis_up$Slope), -rank(abs(cis_down$Slope)))
cis_genes$logrank<-c(log10(rank(cis_up$Slope))+(1e-10), (log10(rank(abs(cis_down$Slope)))*-1)-(1e-10) )
colnames(cis_genes)<-c('Genes', 'Slope', 'Rank', 'LogRank')
print(cis_genes)

cis_genes<-cis_genes[rev(order(cis_genes$Rank)),]

head(cis_genes[rev(order(cis_genes$Rank)),])
tail(cis_genes[rev(order(cis_genes$Rank)),])

addWorksheet(wb1,"CIS_Global")
writeData(wb1, "CIS_Global", cis_genes)

cis_up1<-cis_up[cis_up$Slope>0,]
cis_up1<-cis_up1[rev(order(cis_up1$Slope)),]
addWorksheet(wb1,"CIS_Global_Up")
writeData(wb1, "CIS_Global_Up", cis_up1)

cis_down1<-cis_down[order(cis_down$Slope),]
addWorksheet(wb1,"CIS_Global_Down")
writeData(wb1, "CIS_Global_Down", cis_down1)

#Ranked list to be used for GSEA
CIS_Global<-cis_genes$LogRank
names(CIS_Global)<-cis_genes$Genes
head(CIS_Global)
tail(CIS_Global)

#### GSEA GO ANALYSIS (CIS GLOBAL) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CIS_Global.GO<-gseGO(CIS_Global,
                     OrgDb="org.Mm.eg.db",
                     keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     ont="ALL"#ont="BP",
)

CIS_Global.GO.Plot<-CIS_Global.GO
gseaplot(CIS_Global.GO.Plot, geneSetID=1, by="all")
saveRDS(CIS_Global.GO.Plot, file="EnrichmentPlots/CIS_Global.GO.GSEA.rds")

if(sign(max(CIS_Global.GO$NES))==1) 
{ CIS_Global.GO<-CIS_Global.GO[rev(order(CIS_Global.GO$NES)),]
} else { CIS_Global.GO<-CIS_Global.GO[order(CIS_Global.GO$NES),] }

head(CIS_Global.GO$Description)
tail(CIS_Global.GO$Description)

addWorksheet(wb2,"CIS_Global")
writeData(wb2, "CIS_Global", CIS_Global.GO)

###### GSEA HALLMARK (CIS GLOBAL) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

CIS_Global.HALL<-GSEA(geneList = CIS_Global, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

CIS_Global.HALL.Plot<-CIS_Global.HALL
gseaplot(CIS_Global.HALL.Plot, geneSetID = 1)
saveRDS(CIS_Global.HALL.Plot, file="EnrichmentPlots/CIS_Global.HALL.GSEA.rds")

if(sign(max(CIS_Global.HALL$NES))==1) 
{ CIS_Global.HALL<-CIS_Global.HALL[rev(order(CIS_Global.HALL$NES)),]
} else { CIS_Global.HALL<-CIS_Global.HALL[order(CIS_Global.HALL$NES),] }

head(CIS_Global.HALL$Description)

addWorksheet(wb3,"CIS_Global")
writeData(wb3, "CIS_Global", CIS_Global.HALL)

###### GSEA KEGG PATHWAY (CIS GLOBAL) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to Entrez IDs
gene_ids <- bitr(names(CIS_Global), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
CIS_Global.KEGG<-CIS_Global
names(CIS_Global.KEGG)<-gene_ids$ENTREZID
head(CIS_Global.KEGG)

# Run KEGG GSEA
CIS_Global.KEGG <- gseKEGG(
  geneList = CIS_Global.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

CIS_Global.KEGG.Plot<-CIS_Global.KEGG
gseaplot(CIS_Global.KEGG.Plot, geneSetID = 1)
saveRDS(CIS_Global.KEGG.Plot, file="EnrichmentPlots/CIS_Global.KEGG.GSEA.rds")

if(sign(max(CIS_Global.KEGG$NES))==1) 
{ CIS_Global.KEGG<-CIS_Global.KEGG[rev(order(CIS_Global.KEGG$NES)),]
} else { CIS_Global.KEGG<-CIS_Global.KEGG[order(CIS_Global.KEGG$NES),] }

head(CIS_Global.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-CIS_Global.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(CIS_Global.KEGG$core_enrichment), function(i) {
      #print(CIS_Global.KEGG$core_enrichment[i])
      #print(symbols_list[[i]])
      CIS_Global.KEGG$core_enrichment[i]<-symbols_list[i]
      #print(CIS_Global.KEGG$core_enrichment[i])
      return(CIS_Global.KEGG$core_enrichment[i])
      })

head(xyz)

CIS_Global.modKEGG<-CIS_Global.KEGG
CIS_Global.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(CIS_Global.modKEGG$core_enrichment)


addWorksheet(wb4,"CIS_Global")
writeData(wb4, "CIS_Global", CIS_Global.modKEGG)


######### LINEAR REGRESSION IN A COMPARTMENT SPECIFIC MANNER ###################
#### CIS CORTEX  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset to just cortex (CIS)

#Selecting 190 random points from the whole CIS cortex cross section
#This is done to avoid over-representation of the cortical area.
set.seed(200)
cis_rand_spots.cortex.0h<-sample(rownames(CIS_cortex.0h), sample_n, replace=FALSE)
cis_rand_spots.cortex.12h<-sample(rownames(CIS_cortex.12h), sample_n, replace=FALSE)
cis_rand_spots.cortex.24h<-sample(rownames(CIS_cortex.24h), sample_n, replace=FALSE)
cis_rand_spots.cortex.48h<-sample(rownames(CIS_cortex.48h), sample_n, replace=FALSE)

#Plot to check:
plot(CIS_0h$pos[cis_rand_spots.cortex.0h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Cortex_0h")
plot(CIS_12h$pos[cis_rand_spots.cortex.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Cortex_12h")
plot(CIS_24h$pos[cis_rand_spots.cortex.24h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Cortex_24h")
plot(CIS_48h$pos[cis_rand_spots.cortex.48h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Cortex_48h")

CIS_cortex_gexp <- data.frame(cbind(
  as.matrix(CIS_0h$mat_notlog[genes.shared,cis_rand_spots.cortex.0h]),
  as.matrix(CIS_12h$mat_notlog[genes.shared,cis_rand_spots.cortex.12h]),
  as.matrix(CIS_24h$mat_notlog[genes.shared,cis_rand_spots.cortex.24h]),
  as.matrix(CIS_48h$mat_notlog[genes.shared,cis_rand_spots.cortex.48h]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(cis_rand_spots.cortex.0h)),
  rep(log(12+1), length(cis_rand_spots.cortex.12h)), ## unit of days
  rep(log(24+1), length(cis_rand_spots.cortex.24h)),
  rep(log(48+1), length(cis_rand_spots.cortex.48h))
))  
rownames(meta) <- colnames(CIS_cortex_gexp)
table(meta)

#des = with(meta, model.matrix(~ time)) 
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(CIS_cortex_gexp, design = des) ## fit limma model against time
# vi <- meta$annot == 'cortex'
# des = with(meta[vi,], model.matrix(~ time))
# fit <- lmFit(CIS_gexp[,vi], design = des)
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(CIS_cortex_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

cis_cortex_fc <- top$logFC
names(cis_cortex_fc) <- rownames(top)
cis_cortex_fc<-data.frame(names(cis_cortex_fc) ,cis_cortex_fc)
colnames(cis_cortex_fc)<-c('Gene','Slope')
head(cis_cortex_fc)

cis_cortex_up<-cis_cortex_fc[cis_cortex_fc$Slope>=0,]
head(cis_cortex_up)

cis_cortex_down<-cis_cortex_fc[cis_cortex_fc$Slope<0,]
head(cis_cortex_down)

#Creating a ranked list of genes based on their slope
cis_cortex_genes<-rbind(cis_cortex_up,cis_cortex_down)
print(cis_cortex_genes)
cis_cortex_genes$rank<-c(rank(cis_cortex_up$Slope), -rank(abs(cis_cortex_down$Slope)))
cis_cortex_genes$logrank<-c(log10(rank(cis_cortex_up$Slope))+(1e-10), (log10(rank(abs(cis_cortex_down$Slope)))*-1)-(1e-10) )
colnames(cis_cortex_genes)<-c('Genes', 'Slope', 'Rank', 'LogRank')
print(cis_cortex_genes)

cis_cortex_genes<-cis_cortex_genes[rev(order(cis_cortex_genes$Rank)),]

head(cis_cortex_genes[rev(order(cis_cortex_genes$Rank)),])
tail(cis_cortex_genes[rev(order(cis_cortex_genes$Rank)),])

addWorksheet(wb1,"CIS_Cortex")
writeData(wb1, "CIS_Cortex", cis_cortex_genes)

cis_cortex_up1<-cis_cortex_up[cis_cortex_up$Slope>0,]
cis_cortex_up1<-cis_cortex_up1[rev(order(cis_cortex_up1$Slope)),]
addWorksheet(wb1,"CIS_Cortex_Up")
writeData(wb1, "CIS_Cortex_Up", cis_cortex_up1)

cis_cortex_down1<-cis_cortex_down[order(cis_cortex_down$Slope),]
addWorksheet(wb1,"CIS_Cortex_Down")
writeData(wb1, "CIS_Cortex_Down", cis_cortex_down1)

#Ranked list to be used for GSEA
CIS_Cortex<-cis_cortex_genes$LogRank
names(CIS_Cortex)<-cis_cortex_genes$Genes
head(CIS_Cortex)
tail(CIS_Cortex)

## GSEA GO ANALYSIS (CIS CORTEX) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CIS_Cortex.GO<-gseGO(CIS_Cortex,
                    OrgDb="org.Mm.eg.db",
                    keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    ont="ALL"#ont="BP",
)

CIS_Cortex.GO.Plot<-CIS_Cortex.GO
gseaplot(CIS_Cortex.GO, geneSetID=1, by="all")
saveRDS(CIS_Cortex.GO.Plot, file="EnrichmentPlots/CIS_Cortex.GO.GSEA.rds")

if(sign(max(CIS_Cortex.GO$NES))==1) 
{ CIS_Cortex.GO<-CIS_Cortex.GO[rev(order(CIS_Cortex.GO$NES)),]
} else { CIS_Cortex.GO<-CIS_Cortex.GO[order(CIS_Cortex.GO$NES),] }

head(CIS_Cortex.GO$Description)
tail(CIS_Cortex.GO$Description)

addWorksheet(wb2,"CIS_Cortex")
writeData(wb2, "CIS_Cortex", CIS_Cortex.GO)

###### GSEA HALLMARK (CIS CORTEX) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

CIS_Cortex.HALL<-GSEA(geneList = CIS_Cortex, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

CIS_Cortex.HALL.Plot<-CIS_Cortex.HALL
gseaplot(CIS_Cortex.HALL, geneSetID = 1)
saveRDS(CIS_Cortex.HALL.Plot, file="EnrichmentPlots/CIS_Cortex.HALL.GSEA.rds")

if(sign(max(CIS_Cortex.HALL$NES))==1) 
{ CIS_Cortex.HALL<-CIS_Cortex.HALL[rev(order(CIS_Cortex.HALL$NES)),]
} else { CIS_Cortex.HALL<-CIS_Cortex.HALL[order(CIS_Cortex.HALL$NES),] }

head(CIS_Cortex.HALL$Description)

addWorksheet(wb3,"CIS_Cortex")
writeData(wb3, "CIS_Cortex", CIS_Cortex.HALL)

###### GSEA KEGG PATHWAY (CIS CORTEX) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to Entrez IDs
gene_ids <- bitr(names(CIS_Cortex), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
CIS_Cortex.KEGG<-CIS_Cortex
names(CIS_Cortex.KEGG)<-gene_ids$ENTREZID
head(CIS_Cortex.KEGG)

# Run KEGG GSEA
CIS_Cortex.KEGG <- gseKEGG(
  geneList = CIS_Cortex.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

CIS_Cortex.KEGG.Plot<-CIS_Cortex.KEGG
gseaplot(CIS_Cortex.KEGG, geneSetID = 1)
saveRDS(CIS_Cortex.KEGG.Plot, file="EnrichmentPlots/CIS_Cortex.KEGG_GSEA.rds")

if(sign(max(CIS_Cortex.KEGG$NES))==1) 
{ CIS_Cortex.KEGG<-CIS_Cortex.KEGG[rev(order(CIS_Cortex.KEGG$NES)),]
} else { CIS_Cortex.KEGG<-CIS_Cortex.KEGG[order(CIS_Cortex.KEGG$NES),] }

head(CIS_Cortex.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-CIS_Cortex.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(CIS_Cortex.KEGG$core_enrichment), function(i) {
  #print(CIS_Cortex.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  CIS_Cortex.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(CIS_Cortex.KEGG$core_enrichment[i])
  return(CIS_Cortex.KEGG$core_enrichment[i])
})

head(xyz)

CIS_Cortex.modKEGG<-CIS_Cortex.KEGG
CIS_Cortex.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(CIS_Cortex.modKEGG$core_enrichment)


addWorksheet(wb4,"CIS_Cortex")
writeData(wb4, "CIS_Cortex", CIS_Cortex.modKEGG)

### CIS INTERFACE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Note: INTERFACE is same as Outer Medulla.

#Selecting 190 random points from the whole CIS interface cross section
set.seed(300)
cis_rand_spots.interface.0h<-sample(rownames(CIS_interface.0h), sample_n, replace=FALSE)
cis_rand_spots.interface.12h<-sample(rownames(CIS_interface.12h), sample_n, replace=FALSE)
cis_rand_spots.interface.24h<-sample(rownames(CIS_interface.24h), sample_n, replace=FALSE)
cis_rand_spots.interface.48h<-sample(rownames(CIS_interface.48h), sample_n, replace=FALSE)

#Plot to check:
plot(CIS_0h$pos[cis_rand_spots.interface.0h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Interface_0h")
plot(CIS_12h$pos[cis_rand_spots.interface.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Interface_12h")
plot(CIS_24h$pos[cis_rand_spots.interface.24h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Interface_24h")
plot(CIS_48h$pos[cis_rand_spots.interface.48h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Interface_48h")

CIS_interface_gexp <- data.frame(cbind(
  as.matrix(CIS_0h$mat_notlog[genes.shared,cis_rand_spots.interface.0h]),
  as.matrix(CIS_12h$mat_notlog[genes.shared,cis_rand_spots.interface.12h]),
  as.matrix(CIS_24h$mat_notlog[genes.shared,cis_rand_spots.interface.24h]),
  as.matrix(CIS_48h$mat_notlog[genes.shared,cis_rand_spots.interface.48h]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(cis_rand_spots.interface.0h)),
  rep(log(12+1), length(cis_rand_spots.interface.12h)), ## unit of days
  rep(log(24+1), length(cis_rand_spots.interface.24h)),
  rep(log(48+1), length(cis_rand_spots.interface.48h))
))  
rownames(meta) <- colnames(CIS_interface_gexp)
table(meta)

#vi <- meta$annot == 'interface'
#des = with(meta[vi,], model.matrix(~ time))
#fit <- lmFit(CIS_gexp[,vi], design = des)
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(CIS_interface_gexp, design = des) ## fit limma model against time
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(CIS_interface_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

cis_interface_fc <- top$logFC
names(cis_interface_fc) <- rownames(top)
cis_interface_fc<-data.frame(names(cis_interface_fc) ,cis_interface_fc)
colnames(cis_interface_fc)<-c('Gene','Slope')
head(cis_interface_fc)

cis_interface_up<-cis_interface_fc[cis_interface_fc$Slope>=0,]
head(cis_interface_up)

cis_interface_down<-cis_interface_fc[cis_interface_fc$Slope<0,]
head(cis_interface_down)

#Creating a ranked list of genes based on their slope
cis_interface_genes<-rbind(cis_interface_up,cis_interface_down)
print(cis_interface_genes)
cis_interface_genes$rank<-c(rank(cis_interface_up$Slope), -rank(abs(cis_interface_down$Slope)))
cis_interface_genes$logrank<-c(log10(rank(cis_interface_up$Slope))+(1e-10), (log10(rank(abs(cis_interface_down$Slope)))*-1)-(1e-10) )
colnames(cis_interface_genes)<-c('Genes', 'Slope', 'Rank', 'LogRank')
print(cis_interface_genes)

cis_interface_genes<-cis_interface_genes[rev(order(cis_interface_genes$Rank)),]

cis_interface_up<-cis_interface_fc[cis_interface_fc$Slope>=0,]
head(cis_interface_up)

cis_interface_down<-cis_interface_fc[cis_interface_fc$Slope<0,]
head(cis_interface_down)

head(cis_interface_genes[rev(order(cis_interface_genes$Rank)),])
tail(cis_interface_genes[rev(order(cis_interface_genes$Rank)),])

addWorksheet(wb1,"CIS_Interface")
writeData(wb1, "CIS_Interface", cis_interface_genes)

cis_interface_up1<-cis_interface_up[cis_interface_up$Slope>0,]
cis_interface_up1<-cis_interface_up1[rev(order(cis_interface_up1$Slope)),]
addWorksheet(wb1,"CIS_Interface_Up")
writeData(wb1, "CIS_Interface_Up", cis_interface_up1)

cis_interface_down1<-cis_interface_down[order(cis_interface_down$Slope),]
addWorksheet(wb1,"CIS_Interface_Down")
writeData(wb1, "CIS_Interface_Down", cis_interface_down1)

#Ranked list to be used for GSEA
CIS_Interface<-cis_interface_genes$LogRank
names(CIS_Interface)<-cis_interface_genes$Genes
head(CIS_Interface)
tail(CIS_Interface)

## GSEA GO ANALYSIS (CIS INTERFACE) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CIS_Interface.GO<-gseGO(CIS_Interface,
                     OrgDb="org.Mm.eg.db",
                     keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     ont="ALL"#ont="BP",
)

CIS_Interface.GO.Plot<-CIS_Interface.GO
gseaplot(CIS_Interface.GO, geneSetID=1, by="all")
saveRDS(CIS_Interface.GO.Plot, file="EnrichmentPlots/CIS_Interface.GO_GSEA.rds")

if(sign(max(CIS_Interface.GO$NES))==1) 
{ CIS_Interface.GO<-CIS_Interface.GO[rev(order(CIS_Interface.GO$NES)),]
} else { CIS_Interface.GO<-CIS_Interface.GO[order(CIS_Interface.GO$NES),] }

head(CIS_Interface.GO$Description)
tail(CIS_Interface.GO$Description)

addWorksheet(wb2,"CIS_Interface")
writeData(wb2, "CIS_Interface", CIS_Interface.GO)

###### GSEA HALLMARK (CIS INTERFACE) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

CIS_Interface.HALL<-GSEA(geneList = CIS_Interface, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

CIS_Interface.HALL.Plot<-CIS_Interface.HALL
gseaplot(CIS_Interface.HALL.Plot, geneSetID = 1)
saveRDS(CIS_Interface.HALL.Plot, file="EnrichmentPlots/CIS_Interface.HALL_GSEA.rds")

if(sign(max(CIS_Interface.HALL$NES))==1) 
{ CIS_Interface.HALL<-CIS_Interface.HALL[rev(order(CIS_Interface.HALL$NES)),]
} else { CIS_Interface.HALL<-CIS_Interface.HALL[order(CIS_Interface.HALL$NES),] }

head(CIS_Interface.HALL$Description)

addWorksheet(wb3,"CIS_Interface")
writeData(wb3, "CIS_Interface", CIS_Interface.HALL)

###### GSEA KEGG PATHWAY (CIS INTERFACE) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to Entrez IDs
gene_ids <- bitr(names(CIS_Interface), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
CIS_Interface.KEGG<-CIS_Interface
names(CIS_Interface.KEGG)<-gene_ids$ENTREZID
head(CIS_Interface.KEGG)

# Run KEGG GSEA
CIS_Interface.KEGG <- gseKEGG(
  geneList = CIS_Interface.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

CIS_Interface.KEGG.Plot<-CIS_Interface.KEGG
gseaplot(CIS_Interface.KEGG.Plot, geneSetID = 1)
saveRDS(CIS_Interface.KEGG.Plot, file="EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")

if(sign(max(CIS_Interface.KEGG$NES))==1) 
{ CIS_Interface.KEGG<-CIS_Interface.KEGG[rev(order(CIS_Interface.KEGG$NES)),]
} else { CIS_Interface.KEGG<-CIS_Interface.KEGG[order(CIS_Interface.KEGG$NES),] }

head(CIS_Interface.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-CIS_Interface.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(CIS_Interface.KEGG$core_enrichment), function(i) {
  #print(CIS_Interface.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  CIS_Interface.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(CIS_Interface.KEGG$core_enrichment[i])
  return(CIS_Interface.KEGG$core_enrichment[i])
})

head(xyz)

CIS_Interface.modKEGG<-CIS_Interface.KEGG
CIS_Interface.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(CIS_Interface.modKEGG$core_enrichment)


addWorksheet(wb4,"CIS_Interface")
writeData(wb4, "CIS_Interface", CIS_Interface.modKEGG)


#### CIS MEDULLA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Note: MEDULLA is same as Inner Medulla.

#Selecting 190 random points from the whole CIS medulla cross section
set.seed(400)
cis_rand_spots.medulla.0h<-sample(rownames(CIS_medulla.0h), sample_n, replace=FALSE)
cis_rand_spots.medulla.12h<-sample(rownames(CIS_medulla.12h), sample_n, replace=FALSE)
cis_rand_spots.medulla.24h<-sample(rownames(CIS_medulla.24h), sample_n, replace=FALSE)
cis_rand_spots.medulla.48h<-sample(rownames(CIS_medulla.48h), sample_n, replace=FALSE)

#Plot to check:
plot(CIS_0h$pos[cis_rand_spots.medulla.0h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Medulla_0h")
plot(CIS_12h$pos[cis_rand_spots.medulla.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Medulla_12h")
plot(CIS_24h$pos[cis_rand_spots.medulla.24h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Medulla_24h")
plot(CIS_48h$pos[cis_rand_spots.medulla.48h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="CIS_Medulla_48h")

CIS_medulla_gexp <- data.frame(cbind(
  as.matrix(CIS_0h$mat_notlog[genes.shared,cis_rand_spots.medulla.0h]),
  as.matrix(CIS_12h$mat_notlog[genes.shared,cis_rand_spots.medulla.12h]),
  as.matrix(CIS_24h$mat_notlog[genes.shared,cis_rand_spots.medulla.24h]),
  as.matrix(CIS_48h$mat_notlog[genes.shared,cis_rand_spots.medulla.48h]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(cis_rand_spots.medulla.0h)),
  rep(log(12+1), length(cis_rand_spots.medulla.12h)), ## unit of days
  rep(log(24+1), length(cis_rand_spots.medulla.24h)),
  rep(log(48+1), length(cis_rand_spots.medulla.48h))
))  
rownames(meta) <- colnames(CIS_medulla_gexp)
table(meta)

# vi <- meta$annot == 'medulla'
# des = with(meta[vi,], model.matrix(~ time))
# fit <- lmFit(CIS_gexp[,vi], design = des)
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(CIS_medulla_gexp, design = des) ## fit limma model against time
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(CIS_medulla_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

cis_medulla_fc <- top$logFC
names(cis_medulla_fc) <- rownames(top)
cis_medulla_fc<-data.frame(names(cis_medulla_fc) ,cis_medulla_fc)
colnames(cis_medulla_fc)<-c('Gene','Slope')
head(cis_medulla_fc)

cis_medulla_up<-cis_medulla_fc[cis_medulla_fc$Slope>=0,]
head(cis_medulla_up)

cis_medulla_down<-cis_medulla_fc[cis_medulla_fc$Slope<0,]
head(cis_medulla_down)

#Creating a ranked list of genes based on their slope
cis_medulla_genes<-rbind(cis_medulla_up,cis_medulla_down)
print(cis_medulla_genes)
cis_medulla_genes$rank<-c(rank(cis_medulla_up$Slope), -rank(abs(cis_medulla_down$Slope)))
cis_medulla_genes$logrank<-c(log10(rank(cis_medulla_up$Slope))+(1e-10), (log10(rank(abs(cis_medulla_down$Slope)))*-1)-(1e-10) )
colnames(cis_medulla_genes)<-c('Genes', 'Slope', 'Rank', 'LogRank')
print(cis_medulla_genes)

cis_medulla_genes<-cis_medulla_genes[rev(order(cis_medulla_genes$Rank)),]

head(cis_medulla_genes[rev(order(cis_medulla_genes$Rank)),])
tail(cis_medulla_genes[rev(order(cis_medulla_genes$Rank)),])

addWorksheet(wb1,"CIS_Medulla")
writeData(wb1, "CIS_Medulla", cis_medulla_genes)

cis_medulla_up1<-cis_medulla_up[cis_medulla_up$Slope>0,]
cis_medulla_up1<-cis_medulla_up1[rev(order(cis_medulla_up1$Slope)),]
addWorksheet(wb1,"CIS_Medulla_Up")
writeData(wb1, "CIS_Medulla_Up", cis_medulla_up1)

cis_medulla_down1<-cis_medulla_down[order(cis_medulla_down$Slope),]
addWorksheet(wb1,"CIS_Medulla_Down")
writeData(wb1, "CIS_Medulla_Down", cis_medulla_down1)

#Ranked list to be used for GSEA
CIS_Medulla<-cis_medulla_genes$LogRank
names(CIS_Medulla)<-cis_medulla_genes$Genes
head(CIS_Medulla)
tail(CIS_Medulla)

#GSEA GO ANALYSIS (CIS MEDULLA) ************************************************
CIS_Medulla.GO<-gseGO(CIS_Medulla,
                     OrgDb="org.Mm.eg.db",
                     keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     ont="ALL"#ont="BP",
)

CIS_Medulla.GO.Plot<-CIS_Medulla.GO
gseaplot(CIS_Medulla.GO, geneSetID=1, by="all")
saveRDS(CIS_Medulla.GO.Plot, file="EnrichmentPlots/CIS_Medulla.GO_GSEA.rds")

if(sign(max(CIS_Medulla.GO$NES))==1) 
{ CIS_Medulla.GO<-CIS_Medulla.GO[rev(order(CIS_Medulla.GO$NES)),]
} else { CIS_Medulla.GO<-CIS_Medulla.GO[order(CIS_Medulla.GO$NES),] }

head(CIS_Medulla.GO$Description)
tail(CIS_Medulla.GO$Description)

addWorksheet(wb2,"CIS_Medulla")
writeData(wb2, "CIS_Medulla", CIS_Medulla.GO)

###### GSEA HALLMARK (CIS MEDULLA) *********************************************
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

CIS_Medulla.HALL<-GSEA(geneList = CIS_Medulla, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

CIS_Medulla.HALL.Plot<-CIS_Medulla.HALL
gseaplot(CIS_Medulla.HALL.Plot, geneSetID = 1)
saveRDS(CIS_Medulla.HALL.Plot, file="EnrichmentPlots/CIS_Medulla.HALL_GSEA.rds")

if(sign(max(CIS_Medulla.HALL$NES))==1) 
{ CIS_Medulla.HALL<-CIS_Medulla.HALL[rev(order(CIS_Medulla.HALL$NES)),]
} else { CIS_Medulla.HALL<-CIS_Medulla.HALL[order(CIS_Medulla.HALL$NES),] }

head(CIS_Medulla.HALL$Description)

addWorksheet(wb3,"CIS_Medulla")
writeData(wb3, "CIS_Medulla", CIS_Medulla.HALL)

###### GSEA KEGG PATHWAY (CIS MEDULLA) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to Entrez IDs
gene_ids <- bitr(names(CIS_Medulla), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
CIS_Medulla.KEGG<-CIS_Medulla
names(CIS_Medulla.KEGG)<-gene_ids$ENTREZID
head(CIS_Medulla.KEGG)

# Run KEGG GSEA
CIS_Medulla.KEGG <- gseKEGG(
  geneList = CIS_Medulla.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

CIS_Medulla.KEGG.Plot<-CIS_Medulla.KEGG
gseaplot(CIS_Medulla.KEGG.Plot, geneSetID = 1)
saveRDS(CIS_Medulla.KEGG.Plot, file="EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")

if(sign(max(CIS_Medulla.KEGG$NES))==1) 
{ CIS_Medulla.KEGG<-CIS_Medulla.KEGG[rev(order(CIS_Medulla.KEGG$NES)),]
} else { CIS_Medulla.KEGG<-CIS_Medulla.KEGG[order(CIS_Medulla.KEGG$NES),] }

head(CIS_Medulla.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-CIS_Medulla.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(CIS_Medulla.KEGG$core_enrichment), function(i) {
  #print(CIS_Medulla.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  CIS_Medulla.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(CIS_Medulla.KEGG$core_enrichment[i])
  return(CIS_Medulla.KEGG$core_enrichment[i])
})

head(xyz)

CIS_Medulla.modKEGG<-CIS_Medulla.KEGG
CIS_Medulla.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(CIS_Medulla.modKEGG$core_enrichment)


addWorksheet(wb4,"CIS_Medulla")
writeData(wb4, "CIS_Medulla", CIS_Medulla.modKEGG)

#Saving the genes and GSEA pathways
saveWorkbook(wb1, "Supplementary_Tables/Limma_cpmCIS_DEGs.xlsx", overwrite = TRUE)
saveWorkbook(wb2, "Supplementary_Tables/Limma_cpmCIS_GSEA_GO_ALL_Pathways.xlsx", overwrite = TRUE)
saveWorkbook(wb3, "Supplementary_Tables/Limma_cpmCIS_GSEA_Hallmark_Pathways.xlsx", overwrite = TRUE)
saveWorkbook(wb4, "Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", overwrite = TRUE)


################################################################################
#### LINEAR REGRESSION ANALYSIS FOR WARM ISCHMEMIA-REPERFUSION (AKI) DATASET#### 

wb5<-createWorkbook() #workbook for DEG genelist (AKI sham - 2weeks)
wb6<-createWorkbook() #workbook for gsea GO pathway (AKI sham - 2weeks)
wb7<-createWorkbook() #workbook for gsea HALLMARK pathway (AKI sham - 2weeks)
wb8<-createWorkbook() #workbook for gsea HALLMARK pathway (AKI sham - 2weeks)

genes.shared <- Reduce(intersect, list(rownames(CIS_0h$mat_notlog), rownames(CIS_12h$mat_notlog),
                                       rownames(CIS_24h$mat_notlog), rownames(CIS_48h$mat_notlog),
                                       rownames(AKI_sham$mat_notlog), rownames(AKI_4h$mat_notlog),
                                       rownames(AKI_12h$mat_notlog), rownames(AKI_2d$mat_notlog)) )
#genes.shared<-rownames(CIS_0h$gexp)

#Selecting n=190 random points from the whole AKI tissue cross section
sample_n<-190 #Setting the no. of random spots to be extracted

set.seed(500)
aki_rand_spots.sham<-sample(colnames(AKI_sham$mat_notlog), sample_n, replace=FALSE)
aki_rand_spots.4h<-sample(colnames(AKI_4h$mat_notlog), sample_n, replace=FALSE)
aki_rand_spots.12h<-sample(colnames(AKI_12h$mat_notlog), sample_n, replace=FALSE)
aki_rand_spots.2d<-sample(colnames(AKI_2d$mat_notlog), sample_n, replace=FALSE)

#Plot to check:
plot(AKI_sham$pos[aki_rand_spots.sham,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Global_sham")
plot(AKI_4h$pos[aki_rand_spots.4h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Global_4h")
plot(AKI_12h$pos[aki_rand_spots.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Global_12h")
plot(AKI_2d$pos[aki_rand_spots.2d,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Global_2d")

AKI2d_global_gexp <- data.frame(cbind(
  as.matrix(AKI_sham$mat_notlog[genes.shared,aki_rand_spots.sham]),
  as.matrix(AKI_4h$mat_notlog[genes.shared,aki_rand_spots.4h]),
  as.matrix(AKI_12h$mat_notlog[genes.shared,aki_rand_spots.12h]),
  as.matrix(AKI_2d$mat_notlog[genes.shared,aki_rand_spots.2d])
))

# AKI2d_gexp <- data.frame(cbind(
#   as.matrix(AKI_sham$mat_notlog[genes.shared,]),
#   as.matrix(AKI_4h$mat_notlog[genes.shared,]),
#   as.matrix(AKI_12h$mat_notlog[genes.shared,]),
#   as.matrix(AKI_2d$mat_notlog[genes.shared,])
# ))

meta <- data.frame(time = c(
  rep(log(0+1), length(aki_rand_spots.sham)),
  rep(log(4+1), length(aki_rand_spots.4h)), ## unit of days
  rep(log(12+1), length(aki_rand_spots.12h)),
  rep(log(48+1), length(aki_rand_spots.2d))
))  
rownames(meta) <- colnames(AKI2d_global_gexp)
table(meta)

# meta <- data.frame(time = c(
#   rep(log(0+1), ncol(AKI_sham$mat_notlog)),
#   rep(log(4+1), ncol(AKI_4h$mat_notlog)), ## unit of days
#   rep(log(12+1), ncol(AKI_12h$mat_notlog)),
#   rep(log(48+1), ncol(AKI_2d$mat_notlog))
# ))  
# rownames(meta) <- colnames(AKI2d_global_gexp)
# table(meta)

#des = with(meta, model.matrix(~ time)) 
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(AKI2d_global_gexp, design = des) ## fit limma model against time
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(AKI2d_global_gexp))
top <- top[top$adj.P.Val < 0.05,] ## restrict to significant results only
head(top)

aki_global_fc <- top$logFC
names(aki_global_fc) <- rownames(top)
aki_global_fc<-data.frame(names(aki_global_fc) ,aki_global_fc)
colnames(aki_global_fc)<-c('Gene','Slope')
head(aki_global_fc)

aki_up<-aki_global_fc[aki_global_fc$Slope>=0,]
head(aki_up)

aki_down<-aki_global_fc[aki_global_fc$Slope<0,]
head(aki_down)

#Creating a ranked list of genes based on their slope
aki_genes<-rbind(aki_up,aki_down)
print(aki_genes)
aki_genes$rank<-c(rank(aki_up$Slope), -rank(abs(aki_down$Slope)))
aki_genes$logrank<-c(log10(rank(aki_up$Slope))+(1e-10), (log10(rank(abs(aki_down$Slope)))*-1)-(1e-10) )
colnames(aki_genes)<-c('Genes', 'Slope', 'Rank', 'LogRank')
print(aki_genes)

aki_genes<-aki_genes[rev(order(aki_genes$Rank)),]

head(aki_genes[rev(order(aki_genes$Rank)),])
tail(aki_genes[rev(order(aki_genes$Rank)),])

addWorksheet(wb5,"AKI_Global")
writeData(wb5, "AKI_Global", aki_genes)

aki_up1<-aki_up[aki_up$Slope>0,]
aki_up1<-aki_up1[rev(order(aki_up1$Slope)),]
addWorksheet(wb5,"AKI_Global_Up")
writeData(wb5, "AKI_Global_Up", aki_up1)

aki_down1<-aki_down[order(aki_down$Slope),]
addWorksheet(wb5,"AKI_Global_Down")
writeData(wb5, "AKI_Global_Down", aki_down1)

#Ranked list to be used for GSEA
AKI_Global<-aki_genes$LogRank
names(AKI_Global)<-aki_genes$Genes
head(AKI_Global)
tail(AKI_Global)

#######GSEA GO ANALYSIS (AKI GLOBAL) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AKI_Global.GO<-gseGO(AKI_Global,
                     OrgDb="org.Mm.eg.db",
                     keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     ont="ALL"#ont="BP",
)

AKI_Global.GO.Plot<-AKI_Global.GO
gseaplot(AKI_Global.GO.Plot, geneSetID=1, by="all")
saveRDS(AKI_Global.GO.Plot, file="EnrichmentPlots/AKI_Global.GO_GSEA.rds")

if(sign(max(AKI_Global.GO$NES))==1) 
{ AKI_Global.GO<-AKI_Global.GO[rev(order(AKI_Global.GO$NES)),]
} else { AKI_Global.GO<-AKI_Global.GO[order(AKI_Global.GO$NES),] }

head(AKI_Global.GO$Description)
tail(AKI_Global.GO$Description)

addWorksheet(wb6,"AKI_Global")
writeData(wb6, "AKI_Global", AKI_Global.GO)

###### GSEA HALLMARK (AKI GLOBAL) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

AKI_Global.HALL<-GSEA(geneList = AKI_Global, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

AKI_Global.HALL.Plot<-AKI_Global.HALL
gseaplot(AKI_Global.HALL.Plot, geneSetID = 1)
saveRDS(AKI_Global.HALL.Plot, file="EnrichmentPlots/AKI_Global.HALL_GSEA.rds")

if(sign(max(AKI_Global.HALL$NES))==1) 
{ AKI_Global.HALL<-AKI_Global.HALL[rev(order(AKI_Global.HALL$NES)),]
} else { AKI_Global.HALL<-AKI_Global.HALL[order(AKI_Global.HALL$NES),] }

head(AKI_Global.HALL$Description)

addWorksheet(wb7,"AKI_Global")
writeData(wb7, "AKI_Global", AKI_Global.HALL)

###### GSEA KEGG PATHWAY (AKI GLOBAL) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to Entrez IDs
gene_ids <- bitr(names(AKI_Global), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
AKI_Global.KEGG<-AKI_Global
names(AKI_Global.KEGG)<-gene_ids$ENTREZID
head(AKI_Global.KEGG)

# Run KEGG GSEA
AKI_Global.KEGG <- gseKEGG(
  geneList = AKI_Global.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

AKI_Global.KEGG.Plot<-AKI_Global.KEGG
gseaplot(AKI_Global.KEGG.Plot, geneSetID = 1)
saveRDS(AKI_Global.KEGG.Plot, file="EnrichmentPlots/AKI_Global.KEGG_GSEA.rds")

if(sign(max(AKI_Global.KEGG$NES))==1) 
{ AKI_Global.KEGG<-AKI_Global.KEGG[rev(order(AKI_Global.KEGG$NES)),]
} else { AKI_Global.KEGG<-AKI_Global.KEGG[order(AKI_Global.KEGG$NES),] }

head(AKI_Global.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-AKI_Global.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(AKI_Global.KEGG$core_enrichment), function(i) {
  #print(AKI_Global.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  AKI_Global.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(AKI_Global.KEGG$core_enrichment[i])
  return(AKI_Global.KEGG$core_enrichment[i])
})

head(xyz)

AKI_Global.modKEGG<-AKI_Global.KEGG
AKI_Global.modKEGG$core_enrichment<-xyz

#head(AKI_Global.KEGG$core_enrichment)
head(AKI_Global.modKEGG$core_enrichment)

addWorksheet(wb8,"AKI_Global")
writeData(wb8, "AKI_Global", AKI_Global.modKEGG)


#### AKI CORTEX (Limma) ********************************************************
## subset to just cortex (AKI)

#Selecting 190 random points from the whole AKI cortex cross section
set.seed(600)
aki_rand_spots.cortex.sham<-sample(rownames(AKI_cortex.sham), sample_n, replace=FALSE)
aki_rand_spots.cortex.4h<-sample(rownames(AKI_cortex.4h), sample_n, replace=FALSE)
aki_rand_spots.cortex.12h<-sample(rownames(AKI_cortex.12h), sample_n, replace=FALSE)
aki_rand_spots.cortex.2d<-sample(rownames(AKI_cortex.2d), sample_n, replace=FALSE)

#Plot to check:
plot(AKI_sham$pos[aki_rand_spots.cortex.sham,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Cortex_sham")
plot(AKI_4h$pos[aki_rand_spots.cortex.4h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Cortex_4h")
plot(AKI_12h$pos[aki_rand_spots.cortex.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Cortex_12h")
plot(AKI_2d$pos[aki_rand_spots.cortex.2d,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Cortex_2d")

AKI_cortex_gexp <- data.frame(cbind(
  as.matrix(AKI_sham$mat_notlog[genes.shared,aki_rand_spots.cortex.sham]),
  as.matrix(AKI_4h$mat_notlog[genes.shared,aki_rand_spots.cortex.4h]),
  as.matrix(AKI_12h$mat_notlog[genes.shared,aki_rand_spots.cortex.12h]),
  as.matrix(AKI_2d$mat_notlog[genes.shared,aki_rand_spots.cortex.2d]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(aki_rand_spots.cortex.sham)),
  rep(log(4+1), length(aki_rand_spots.cortex.4h)), ## unit of days
  rep(log(12+1), length(aki_rand_spots.cortex.12h)),
  rep(log(48+1), length(aki_rand_spots.cortex.2d))
))  
rownames(meta) <- colnames(AKI_cortex_gexp)
table(meta)

#des = with(meta, model.matrix(~ time)) 
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(AKI_cortex_gexp, design = des) ## fit limma model against time
# vi <- meta$annot == 'cortex'
# des = with(meta[vi,], model.matrix(~ time))
# fit <- lmFit(CIS_gexp[,vi], design = des)
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(AKI_cortex_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

aki_cortex_fc <- top$logFC
names(aki_cortex_fc) <- rownames(top)
aki_cortex_fc<-data.frame(names(aki_cortex_fc) ,aki_cortex_fc)
colnames(aki_cortex_fc)<-c('Gene','Slope')
head(aki_cortex_fc)

aki_cortex_up<-aki_cortex_fc[aki_cortex_fc$Slope>=0,]
head(aki_cortex_up)

aki_cortex_down<-aki_cortex_fc[aki_cortex_fc$Slope<0,]
head(aki_cortex_down)

#Creating a ranked list of genes based on their slope
aki_cortex_genes<-rbind(aki_cortex_up,aki_cortex_down)
print(aki_cortex_genes)
aki_cortex_genes$rank<-c(rank(aki_cortex_up$Slope), -rank(abs(aki_cortex_down$Slope)))
aki_cortex_genes$logrank<-c(log10(rank(aki_cortex_up$Slope))+(1e-10), (log10(rank(abs(aki_cortex_down$Slope)))*-1)-(1e-10) )
colnames(aki_cortex_genes)<-c('Genes', 'Slope', 'Rank', 'LogRank')
print(aki_cortex_genes)

aki_cortex_genes<-aki_cortex_genes[rev(order(aki_cortex_genes$Rank)),]

head(aki_cortex_genes[rev(order(aki_cortex_genes$Rank)),])
tail(aki_cortex_genes[rev(order(aki_cortex_genes$Rank)),])

addWorksheet(wb5,"AKI_Cortex")
writeData(wb5, "AKI_Cortex", aki_cortex_genes)

aki_cortex_up1<-aki_cortex_up[aki_cortex_up$Slope>0,]
aki_cortex_up1<-aki_cortex_up1[rev(order(aki_cortex_up1$Slope)),]
addWorksheet(wb5,"AKI_Cortex_Up")
writeData(wb5, "AKI_Cortex_Up", aki_cortex_up1)

aki_cortex_down1<-aki_cortex_down[order(aki_cortex_down$Slope),]
addWorksheet(wb5,"AKI_Cortex_Down")
writeData(wb5, "AKI_Cortex_Down", aki_cortex_down1)

#Ranked list to be used for GSEA
AKI_Cortex<-aki_cortex_genes$LogRank
names(AKI_Cortex)<-aki_cortex_genes$Genes
head(AKI_Cortex)
tail(AKI_Cortex)

## GSEA GO ANALYSIS (AKI CORTEX) ***********************************************
AKI_Cortex.GO<-gseGO(AKI_Cortex,
                     OrgDb="org.Mm.eg.db",
                     keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     ont="ALL"#ont="BP",
)

AKI_Cortex.GO.Plot<-AKI_Cortex.GO
gseaplot(AKI_Cortex.GO, geneSetID=1, by="all")
saveRDS(AKI_Cortex.GO.Plot, file="EnrichmentPlots/AKI_Cortex.GO_GSEA.rds")

if(sign(max(AKI_Cortex.GO$NES))==1) 
{ AKI_Cortex.GO<-AKI_Cortex.GO[rev(order(AKI_Cortex.GO$NES)),]
} else { AKI_Cortex.GO<-AKI_Cortex.GO[order(AKI_Cortex.GO$NES),] }

head(AKI_Cortex.GO$Description)
tail(AKI_Cortex.GO$Description)

addWorksheet(wb6,"AKI_Cortex")
writeData(wb6, "AKI_Cortex", AKI_Cortex.GO)

###### GSEA HALLMARK (AKI CORTEX) **********************************************
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

AKI_Cortex.HALL<-GSEA(geneList = AKI_Cortex, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

AKI_Cortex.HALL.Plot<-AKI_Cortex.HALL
gseaplot(AKI_Cortex.HALL, geneSetID = 1)
saveRDS(AKI_Cortex.HALL.Plot, file="EnrichmentPlots/AKI_Cortex.HALL_GSEA.rds")

if(sign(max(AKI_Cortex.HALL$NES))==1) 
{ AKI_Cortex.HALL<-AKI_Cortex.HALL[rev(order(AKI_Cortex.HALL$NES)),]
} else { AKI_Cortex.HALL<-AKI_Cortex.HALL[order(AKI_Cortex.HALL$NES),] }

head(AKI_Cortex.HALL$Description)

addWorksheet(wb7,"AKI_Cortex")
writeData(wb7, "AKI_Cortex", AKI_Cortex.HALL)

###### GSEA KEGG PATHWAY (AKI CORTEX) ******************************************
# Convert to Entrez IDs
gene_ids <- bitr(names(AKI_Cortex), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
AKI_Cortex.KEGG<-AKI_Cortex
names(AKI_Cortex.KEGG)<-gene_ids$ENTREZID
head(AKI_Cortex.KEGG)

# Run KEGG GSEA
AKI_Cortex.KEGG <- gseKEGG(
  geneList = AKI_Cortex.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

AKI_Cortex.KEGG.Plot<-AKI_Cortex.KEGG
gseaplot(AKI_Cortex.KEGG, geneSetID = 1)
saveRDS(AKI_Cortex.KEGG.Plot, file="EnrichmentPlots/AKI_Cortex.KEGG_GSEA.rds")

if(sign(max(AKI_Cortex.KEGG$NES))==1) 
{ AKI_Cortex.KEGG<-AKI_Cortex.KEGG[rev(order(AKI_Cortex.KEGG$NES)),]
} else { AKI_Cortex.KEGG<-AKI_Cortex.KEGG[order(AKI_Cortex.KEGG$NES),] }

head(AKI_Cortex.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-AKI_Cortex.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(AKI_Cortex.KEGG$core_enrichment), function(i) {
  #print(AKI_Cortex.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  AKI_Cortex.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(AKI_Cortex.KEGG$core_enrichment[i])
  return(AKI_Cortex.KEGG$core_enrichment[i])
})

head(xyz)

AKI_Cortex.modKEGG<-AKI_Cortex.KEGG
AKI_Cortex.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(AKI_Cortex.modKEGG$core_enrichment)


addWorksheet(wb8,"AKI_Cortex")
writeData(wb8, "AKI_Cortex", AKI_Cortex.modKEGG)

#### AKI INTERFACE (Limma) ********************************************************
## subset to just interface (AKI)

#Selecting 190 random points from the whole AKI interface cross section
set.seed(700)
aki_rand_spots.interface.sham<-sample(rownames(AKI_interface.sham), sample_n, replace=FALSE)
aki_rand_spots.interface.4h<-sample(rownames(AKI_interface.4h), sample_n, replace=FALSE)
aki_rand_spots.interface.12h<-sample(rownames(AKI_interface.12h), sample_n, replace=FALSE)
aki_rand_spots.interface.2d<-sample(rownames(AKI_interface.2d), sample_n, replace=FALSE)

#Plot to check:
plot(AKI_sham$pos[aki_rand_spots.interface.sham,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Interface_sham")
plot(AKI_4h$pos[aki_rand_spots.interface.4h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Interface_4h")
plot(AKI_12h$pos[aki_rand_spots.interface.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Interface_12h")
plot(AKI_2d$pos[aki_rand_spots.interface.2d,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Interface_2d")

AKI_interface_gexp <- data.frame(cbind(
  as.matrix(AKI_sham$mat_notlog[genes.shared,aki_rand_spots.interface.sham]),
  as.matrix(AKI_4h$mat_notlog[genes.shared,aki_rand_spots.interface.4h]),
  as.matrix(AKI_12h$mat_notlog[genes.shared,aki_rand_spots.interface.12h]),
  as.matrix(AKI_2d$mat_notlog[genes.shared,aki_rand_spots.interface.2d]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(aki_rand_spots.interface.sham)),
  rep(log(4+1), length(aki_rand_spots.interface.4h)), ## unit of days
  rep(log(12+1), length(aki_rand_spots.interface.12h)),
  rep(log(48+1), length(aki_rand_spots.interface.2d))
))  
rownames(meta) <- colnames(AKI_interface_gexp)
table(meta)

#des = with(meta, model.matrix(~ time)) 
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(AKI_interface_gexp, design = des) ## fit limma model against time
# vi <- meta$annot == 'interface'
# des = with(meta[vi,], model.matrix(~ time))
# fit <- lmFit(CIS_gexp[,vi], design = des)
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(AKI_interface_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

aki_interface_fc <- top$logFC
names(aki_interface_fc) <- rownames(top)
aki_interface_fc<-data.frame(names(aki_interface_fc) ,aki_interface_fc)
colnames(aki_interface_fc)<-c('Gene','Slope')
head(aki_interface_fc)

aki_interface_up<-aki_interface_fc[aki_interface_fc$Slope>=0,]
head(aki_interface_up)

aki_interface_down<-aki_interface_fc[aki_interface_fc$Slope<0,]
head(aki_interface_down)

#Creating a ranked list of genes based on their slope
aki_interface_genes<-rbind(aki_interface_up,aki_interface_down)
print(aki_interface_genes)
aki_interface_genes$rank<-c(rank(aki_interface_up$Slope), -rank(abs(aki_interface_down$Slope)))
aki_interface_genes$logrank<-c(log10(rank(aki_interface_up$Slope))+(1e-10), (log10(rank(abs(aki_interface_down$Slope)))*-1)-(1e-10) )
colnames(aki_interface_genes)<-c('Genes', 'Slope', 'Rank', 'LogRank')
print(aki_interface_genes)

aki_interface_genes<-aki_interface_genes[rev(order(aki_interface_genes$Rank)),]

head(aki_interface_genes[rev(order(aki_interface_genes$Rank)),])
tail(aki_interface_genes[rev(order(aki_interface_genes$Rank)),])

addWorksheet(wb5,"AKI_Interface")
writeData(wb5, "AKI_Interface", aki_interface_genes)

aki_interface_up1<-aki_interface_up[aki_interface_up$Slope>0,]
aki_interface_up1<-aki_interface_up1[rev(order(aki_interface_up1$Slope)),]
addWorksheet(wb5,"AKI_Interface_Up")
writeData(wb5, "AKI_Interface_Up", aki_interface_up1)

aki_interface_down1<-aki_interface_down[order(aki_interface_down$Slope),]
addWorksheet(wb5,"AKI_Interface_Down")
writeData(wb5, "AKI_Interface_Down", aki_interface_down1)

#Ranked list to be used for GSEA
AKI_Interface<-aki_interface_genes$LogRank
names(AKI_Interface)<-aki_interface_genes$Genes
head(AKI_Interface)
tail(AKI_Interface)

## GSEA GO ANALYSIS (AKI INTERFACE) ***********************************************
AKI_Interface.GO<-gseGO(AKI_Interface,
                     OrgDb="org.Mm.eg.db",
                     keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     ont="ALL"#ont="BP",
)

AKI_Interface.GO.Plot<-AKI_Interface.GO
gseaplot(AKI_Interface.GO, geneSetID=1, by="all")
saveRDS(AKI_Interface.GO.Plot, file="EnrichmentPlots/AKI_Interface.GO_GSEA.rds")

if(sign(max(AKI_Interface.GO$NES))==1) 
{ AKI_Interface.GO<-AKI_Interface.GO[rev(order(AKI_Interface.GO$NES)),]
} else { AKI_Interface.GO<-AKI_Interface.GO[order(AKI_Interface.GO$NES),] }

head(AKI_Interface.GO$Description)
tail(AKI_Interface.GO$Description)

addWorksheet(wb6,"AKI_Interface")
writeData(wb6, "AKI_Interface", AKI_Interface.GO)

###### GSEA HALLMARK (AKI INTERFACE) **********************************************
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

AKI_Interface.HALL<-GSEA(geneList = AKI_Interface, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

AKI_Interface.HALL.Plot<-AKI_Interface.HALL
gseaplot(AKI_Interface.HALL, geneSetID = 1)
saveRDS(AKI_Interface.HALL.Plot, file="EnrichmentPlots/AKI_Interface.HALL_GSEA.rds")

if(sign(max(AKI_Interface.HALL$NES))==1) 
{ AKI_Interface.HALL<-AKI_Interface.HALL[rev(order(AKI_Interface.HALL$NES)),]
} else { AKI_Interface.HALL<-AKI_Interface.HALL[order(AKI_Interface.HALL$NES),] }

head(AKI_Interface.HALL$Description)

addWorksheet(wb7,"AKI_Interface")
writeData(wb7, "AKI_Interface", AKI_Interface.HALL)

###### GSEA KEGG PATHWAY (AKI INTERFACE) ******************************************
# Convert to Entrez IDs
gene_ids <- bitr(names(AKI_Interface), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
AKI_Interface.KEGG<-AKI_Interface
names(AKI_Interface.KEGG)<-gene_ids$ENTREZID
head(AKI_Interface.KEGG)

# Run KEGG GSEA
AKI_Interface.KEGG <- gseKEGG(
  geneList = AKI_Interface.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

AKI_Interface.KEGG.Plot<-AKI_Interface.KEGG
gseaplot(AKI_Interface.KEGG, geneSetID = 1)
saveRDS(AKI_Interface.KEGG.Plot, file="EnrichmentPlots/AKI_Interface.KEGG_GSEA.rds")

if(sign(max(AKI_Interface.KEGG$NES))==1) 
{ AKI_Interface.KEGG<-AKI_Interface.KEGG[rev(order(AKI_Interface.KEGG$NES)),]
} else { AKI_Interface.KEGG<-AKI_Interface.KEGG[order(AKI_Interface.KEGG$NES),] }

head(AKI_Interface.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-AKI_Interface.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(AKI_Interface.KEGG$core_enrichment), function(i) {
  #print(AKI_Interface.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  AKI_Interface.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(AKI_Interface.KEGG$core_enrichment[i])
  return(AKI_Interface.KEGG$core_enrichment[i])
})

head(xyz)

AKI_Interface.modKEGG<-AKI_Interface.KEGG
AKI_Interface.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(AKI_Interface.modKEGG$core_enrichment)


addWorksheet(wb8,"AKI_Interface")
writeData(wb8, "AKI_Interface", AKI_Interface.modKEGG)

#### AKI MEDULLA (Limma) ********************************************************
## subset to just interface (AKI)

#Selecting 190 random points from the whole AKI medulla cross section
set.seed(800)
aki_rand_spots.medulla.sham<-sample(rownames(AKI_medulla.sham), sample_n, replace=FALSE)
aki_rand_spots.medulla.4h<-sample(rownames(AKI_medulla.4h), sample_n, replace=FALSE)
aki_rand_spots.medulla.12h<-sample(rownames(AKI_medulla.12h), sample_n, replace=FALSE)
aki_rand_spots.medulla.2d<-sample(rownames(AKI_medulla.2d), sample_n, replace=FALSE)

#Plot to check:
plot(AKI_sham$pos[aki_rand_spots.medulla.sham,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Medulla_sham")
plot(AKI_4h$pos[aki_rand_spots.medulla.4h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Medulla_4h")
plot(AKI_12h$pos[aki_rand_spots.medulla.12h,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Medulla_12h")
plot(AKI_2d$pos[aki_rand_spots.medulla.2d,], pch = 21, bg = "black", col = "black", cex = 0.5, main="AKI_Medulla_2d")

AKI_medulla_gexp <- data.frame(cbind(
  as.matrix(AKI_sham$mat_notlog[genes.shared,aki_rand_spots.medulla.sham]),
  as.matrix(AKI_4h$mat_notlog[genes.shared,aki_rand_spots.medulla.4h]),
  as.matrix(AKI_12h$mat_notlog[genes.shared,aki_rand_spots.medulla.12h]),
  as.matrix(AKI_2d$mat_notlog[genes.shared,aki_rand_spots.medulla.2d]) ))

meta <- data.frame(time = c(
  rep(log(0+1), length(aki_rand_spots.medulla.sham)),
  rep(log(4+1), length(aki_rand_spots.medulla.4h)), ## unit of days
  rep(log(12+1), length(aki_rand_spots.medulla.12h)),
  rep(log(48+1), length(aki_rand_spots.medulla.2d))
))  
rownames(meta) <- colnames(AKI_medulla_gexp)
table(meta)

#des = with(meta, model.matrix(~ time)) 
des = model.matrix(~ time, data = meta) #creating the design matrix
fit <- lmFit(AKI_medulla_gexp, design = des) ## fit limma model against time
# vi <- meta$annot == 'medulla'
# des = with(meta[vi,], model.matrix(~ time))
# fit <- lmFit(CIS_gexp[,vi], design = des)
head(fit$coefficients)
fit <- eBayes(fit)
top <- topTable(fit, number=nrow(AKI_medulla_gexp))
top <- top[top$adj.P.Val < 0.05,]
head(top)

aki_medulla_fc <- top$logFC
names(aki_medulla_fc) <- rownames(top)
aki_medulla_fc<-data.frame(names(aki_medulla_fc) ,aki_medulla_fc)
colnames(aki_medulla_fc)<-c('Gene','Slope')
head(aki_medulla_fc)

aki_medulla_up<-aki_medulla_fc[aki_medulla_fc$Slope>=0,]
head(aki_medulla_up)

aki_medulla_down<-aki_medulla_fc[aki_medulla_fc$Slope<0,]
head(aki_medulla_down)

#Creating a ranked list of genes based on their slope
aki_medulla_genes<-rbind(aki_medulla_up,aki_medulla_down)
print(aki_medulla_genes)
aki_medulla_genes$rank<-c(rank(aki_medulla_up$Slope), -rank(abs(aki_medulla_down$Slope)))
aki_medulla_genes$logrank<-c(log10(rank(aki_medulla_up$Slope))+(1e-10), (log10(rank(abs(aki_medulla_down$Slope)))*-1)-(1e-10) )
colnames(aki_medulla_genes)<-c('Genes', 'Slope', 'Rank', 'LogRank')
print(aki_medulla_genes)

aki_medulla_genes<-aki_medulla_genes[rev(order(aki_medulla_genes$Rank)),]

head(aki_medulla_genes[rev(order(aki_medulla_genes$Rank)),])
tail(aki_medulla_genes[rev(order(aki_medulla_genes$Rank)),])

addWorksheet(wb5,"AKI_Medulla")
writeData(wb5, "AKI_Medulla", aki_medulla_genes)

aki_medulla_up1<-aki_medulla_up[aki_medulla_up$Slope>0,]
aki_medulla_up1<-aki_medulla_up1[rev(order(aki_medulla_up1$Slope)),]
addWorksheet(wb5,"AKI_Medulla_Up")
writeData(wb5, "AKI_Medulla_Up", aki_medulla_up1)

aki_medulla_down1<-aki_medulla_down[order(aki_medulla_down$Slope),]
addWorksheet(wb5,"AKI_Medulla_Down")
writeData(wb5, "AKI_Medulla_Down", aki_medulla_down1)

#Ranked list to be used for GSEA
AKI_Medulla<-aki_medulla_genes$LogRank
names(AKI_Medulla)<-aki_medulla_genes$Genes
head(AKI_Medulla)
tail(AKI_Medulla)

## GSEA GO ANALYSIS (AKI MEDULLA) ***********************************************
AKI_Medulla.GO<-gseGO(AKI_Medulla,
                        OrgDb="org.Mm.eg.db",
                        keyType="SYMBOL", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                        ont="ALL"#ont="BP",
)

AKI_Medulla.GO.Plot<-AKI_Medulla.GO
gseaplot(AKI_Medulla.GO, geneSetID=1, by="all")
saveRDS(AKI_Medulla.GO.Plot, file="EnrichmentPlots/AKI_Medulla.GO_GSEA.rds")

if(sign(max(AKI_Medulla.GO$NES))==1) 
{ AKI_Medulla.GO<-AKI_Medulla.GO[rev(order(AKI_Medulla.GO$NES)),]
} else { AKI_Medulla.GO<-AKI_Medulla.GO[order(AKI_Medulla.GO$NES),] }

head(AKI_Medulla.GO$Description)
tail(AKI_Medulla.GO$Description)

addWorksheet(wb6,"AKI_Medulla")
writeData(wb6, "AKI_Medulla", AKI_Medulla.GO)

###### GSEA HALLMARK (AKI MEDULLA) **********************************************
#Creating background gene list for GSEA 
gene_sets<-msigdbr(species = "mouse", category = "H")
background <- gene_sets[, c("gs_name", "gene_symbol")]

AKI_Medulla.HALL<-GSEA(geneList = AKI_Medulla, TERM2GENE = background, pvalueCutoff = 0.05, pAdjustMethod = "BH",)

AKI_Medulla.HALL.Plot<-AKI_Medulla.HALL
gseaplot(AKI_Medulla.HALL, geneSetID = 1)
saveRDS(AKI_Medulla.HALL.Plot, file="EnrichmentPlots/AKI_Medulla.HALL_GSEA.rds")

if(sign(max(AKI_Medulla.HALL$NES))==1) 
{ AKI_Medulla.HALL<-AKI_Medulla.HALL[rev(order(AKI_Medulla.HALL$NES)),]
} else { AKI_Medulla.HALL<-AKI_Medulla.HALL[order(AKI_Medulla.HALL$NES),] }

head(AKI_Medulla.HALL$Description)

addWorksheet(wb7,"AKI_Medulla")
writeData(wb7, "AKI_Medulla", AKI_Medulla.HALL)

###### GSEA KEGG PATHWAY (AKI MEDULLA) ******************************************
# Convert to Entrez IDs
gene_ids <- bitr(names(AKI_Medulla), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=FALSE)
AKI_Medulla.KEGG<-AKI_Medulla
names(AKI_Medulla.KEGG)<-gene_ids$ENTREZID
head(AKI_Medulla.KEGG)

# Run KEGG GSEA
AKI_Medulla.KEGG <- gseKEGG(
  geneList = AKI_Medulla.KEGG,          # Ranked gene list with Entrez IDs
  organism = "mmu",                 # "hsa" for human; use "mmu" for mouse, etc.
  pvalueCutoff = 0.05,              # P-value threshold
  verbose = FALSE                   # Suppress verbose output
)

AKI_Medulla.KEGG.Plot<-AKI_Medulla.KEGG
gseaplot(AKI_Medulla.KEGG, geneSetID = 1)
saveRDS(AKI_Medulla.KEGG.Plot, file="EnrichmentPlots/AKI_Medulla.KEGG_GSEA.rds")

if(sign(max(AKI_Medulla.KEGG$NES))==1) 
{ AKI_Medulla.KEGG<-AKI_Medulla.KEGG[rev(order(AKI_Medulla.KEGG$NES)),]
} else { AKI_Medulla.KEGG<-AKI_Medulla.KEGG[order(AKI_Medulla.KEGG$NES),] }

head(AKI_Medulla.KEGG$Description)

#To convert the entrezID back to gene symbols
entrez_ids<-AKI_Medulla.KEGG$core_enrichment

split_ids <- strsplit(entrez_ids, "/") #To separate the entrezIDs 

gene_symbol<-lapply(split_ids, function(entrez_sublist) {
  bitr(entrez_sublist, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
})

symbols_list <- lapply(gene_symbol, function(df) df$SYMBOL)
head(symbols_list)

xyz<-lapply(seq_along(AKI_Medulla.KEGG$core_enrichment), function(i) {
  #print(AKI_Medulla.KEGG$core_enrichment[i])
  #print(symbols_list[[i]])
  AKI_Medulla.KEGG$core_enrichment[i]<-symbols_list[i]
  #print(AKI_Medulla.KEGG$core_enrichment[i])
  return(AKI_Medulla.KEGG$core_enrichment[i])
})

head(xyz)

AKI_Medulla.modKEGG<-AKI_Medulla.KEGG
AKI_Medulla.modKEGG$core_enrichment<-xyz

#head(CIS_Global.KEGG$core_enrichment)
head(AKI_Medulla.modKEGG$core_enrichment)

addWorksheet(wb8,"AKI_Medulla")
writeData(wb8, "AKI_Medulla", AKI_Medulla.modKEGG)


#Saving the genes and GSEA pathways
saveWorkbook(wb5, "Supplementary_Tables/Limma_cpmAKI2d_DEGs.xlsx", overwrite = TRUE)
saveWorkbook(wb6, "Supplementary_Tables/Limma_cpmAKI2d_GSEA_GO_ALL_Pathways.xlsx", overwrite = TRUE)
saveWorkbook(wb7, "Supplementary_Tables/Limma_cpmAKI2d_GSEA_HALLMARK_Pathways.xlsx", overwrite = TRUE)
saveWorkbook(wb8, "Supplementary_Tables/Limma_cpmAKI2d_GSEA_KEGG_Pathways.xlsx", overwrite = TRUE)

