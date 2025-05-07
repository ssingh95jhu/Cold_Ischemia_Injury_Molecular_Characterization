#Clean environment
rm(list=ls(all.names=TRUE))
gc()

load("data/CIS_data.RData")

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(cowplot)
library(viridis) #color blind friendly
library(RColorBrewer) #for RColorbrewer

#CPM Normalization
CIS_0h_cd<-MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_12h_cd<-MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_24h_cd<-MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_48h_cd<-MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
  
#Extracting information
CIS_Medulla<-read.xlsx("Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Medulla")
CIS_Medulla$Description<-gsub(" - Mus.*","",CIS_Medulla$Description)

#### PRINCIPAL COMPONENT PLOTS #################################################
## ASSUMING PCs ARE NOT SAME ACROSS COLD ISCHEMIA TIME POINTS
##~~~~~~~~~~~~~~~~~~~~~~~ OXPHOS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#CIS_OXPHOS_Med_genes<-strsplit(CIS_Medulla[grep("Oxidative phosphorylation", CIS_Medulla$Description ),]$core_enrichment,", ")[[1]][1:10]
CIS_OXPHOS_Med_genes<-strsplit(CIS_Medulla[grep("Oxidative phosphorylation", CIS_Medulla$Description ),]$core_enrichment,", ")[[1]]

##CIS 0hrs
CIS_0h_Med_OXPHOS<-CIS_0h_cd[CIS_OXPHOS_Med_genes,]

set.seed(100)
pcs <- MUDAN::getPcs(CIS_0h_Med_OXPHOS,
                     nGenes=length(CIS_0h_Med_OXPHOS), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X, Y=CIS_0h$pos[names(PC1),]$Y,PC1=PC1)
ggplot(df1, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 0h_OXPHOS")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z1<-PC1
Z1_max<-max(PC1)
Z1_min<-min(PC1)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X, Y=CIS_0h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df1, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 0h_OXPHOS")

##CIS 12hrs
CIS_12h_Med_OXPHOS<-CIS_12h_cd[CIS_OXPHOS_Med_genes,]
set.seed(100)
pcs <- MUDAN::getPcs(CIS_12h_Med_OXPHOS,
                     nGenes=length(CIS_12h_Med_OXPHOS), 
                     nPcs=5,
                     verbose=FALSE) 

#Normalize PC1
PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X, Y=CIS_12h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df2, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 12h_OXPHOS")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z2<-PC1
Z2_max<-max(PC1)
Z2_min<-min(PC1)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X, Y=CIS_12h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df2, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 12h_OXPHOS")

##CIS 24hrs
CIS_24h_Med_OXPHOS<-CIS_24h_cd[CIS_OXPHOS_Med_genes,]
set.seed(300)
pcs <- MUDAN::getPcs(CIS_24h_Med_OXPHOS,
                     nGenes=length(CIS_24h_Med_OXPHOS), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X, Y=CIS_24h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df3, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 24h_OXPHOS")

#zscore of PC!
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z3<-PC1
Z3_max<-max(PC1)
Z3_min<-min(PC1)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X, Y=CIS_24h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df3, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 24h_OXPHOS")

##CIS_48hrs
CIS_48h_Med_OXPHOS<-CIS_48h_cd[CIS_OXPHOS_Med_genes,]
set.seed(500)
pcs <- MUDAN::getPcs(CIS_48h_Med_OXPHOS,
                     nGenes=length(CIS_48h_Med_OXPHOS), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X, Y=CIS_48h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df4, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 48h_OXPHOS")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z4<-PC1
Z4_max<-max(PC1)
Z4_min<-min(PC1)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X, Y=CIS_48h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df4, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 48h_OXPHOS")

#Plot with common scale:
Zmax<-quantile(c(Z1,Z2,Z3,Z4), 0.99)
Zmin<-min(Z1_min, Z2_min, Z3_min, Z4_min)

color_scale <- scale_color_gradientn(colors=viridis(5), values = c(0, 0.8, 1),
                                     limits = c(Zmin, Zmax), oob = scales::squish)

ggplot(df1, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 0h_OXPHOS") + color_scale
ggplot(df2, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 12h_OXPHOS") + color_scale
ggplot(df3, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 24h_OXPHOS") + color_scale
ggplot(df4, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 48h_OXPHOS") + color_scale

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~ GLYCOLYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Extracting information
CIS_Medulla<-read.xlsx("Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Medulla")
CIS_Medulla$Description<-gsub(" - Mus.*","",CIS_Medulla$Description)


CIS_Glycolysis_Med_genes<-strsplit(CIS_Medulla[grep("Glycolysis / Gluconeogenesis", CIS_Medulla$Description ),]$core_enrichment,", ")[[1]]

set.seed(100)
##CIS 0hrs
CIS_0h_Med_Glycolysis<-CIS_0h_cd[CIS_Glycolysis_Med_genes,]

pcs <- MUDAN::getPcs(CIS_0h_Med_Glycolysis,
                     nGenes=length(CIS_0h_Med_Glycolysis), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X,
                Y=CIS_0h$pos[names(PC1),]$Y,
                PC1=PC1)

ggplot(df1, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 0h_Glycolysis")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z1<-PC1
Z1_max<-max(PC1)
Z1_min<-min(PC1)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X, Y=CIS_0h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df1, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 0h_Glycolysis")

set.seed(200)
##CIS 12hrs
CIS_12h_Med_Glycolysis<-CIS_12h_cd[CIS_Glycolysis_Med_genes,]

pcs <- MUDAN::getPcs(CIS_12h_Med_Glycolysis,
                     nGenes=length(CIS_12h_Med_Glycolysis), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X, Y=CIS_12h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df2, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 12h_Glycolysis")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z2<-PC1
Z2_max<-max(PC1)
Z2_min<-min(PC1)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X, Y=CIS_12h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df2, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 12h_Glycolysis")

##CIS 24hrs
CIS_24h_Med_Glycolysis<-CIS_24h_cd[CIS_Glycolysis_Med_genes,]
set.seed(300)
pcs <- MUDAN::getPcs(CIS_24h_Med_Glycolysis,
                     nGenes=length(CIS_24h_Med_Glycolysis), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X, Y=CIS_24h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df3, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 24h_Glycolysis")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z3<-PC1
Z3_max<-max(PC1)
Z3_min<-min(PC1)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X, Y=CIS_24h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df3, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 24h_Glycolysis")

##CIS 48hrs
CIS_48h_Med_Glycolysis<-CIS_48h_cd[CIS_Glycolysis_Med_genes,]
set.seed(500)
pcs <- MUDAN::getPcs(CIS_48h_Med_Glycolysis,
                     nGenes=length(CIS_48h_Med_Glycolysis), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X, Y=CIS_48h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df4, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 48h_Glycolysis")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z4<-PC1
Z4_max<-max(PC1)
Z4_min<-min(PC1)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X, Y=CIS_48h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df4, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 48h_Glycolysis")

#Plot with common scale:
Zmax<-quantile(c(Z1,Z2,Z3,Z4), 0.99)
Zmin<-min(Z1_min, Z2_min, Z3_min, Z4_min)
Zmin

color_scale <- scale_color_gradientn(colors=viridis(5), values = c(0, 0.8, 1),
                                     limits = c(Zmin, Zmax), oob = scales::squish)

ggplot(df1, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 0h_Glycolysis") + color_scale
ggplot(df2, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 12h_Glycolysis") + color_scale
ggplot(df3, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 24h_Glycolysis") + color_scale
ggplot(df4, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 48h_Glycolysis") + color_scale


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ROS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CIS_ROS_Med_genes<-strsplit(CIS_Medulla[grep("Chemical carcinogenesis - reactive oxygen species", CIS_Medulla$Description ),]$core_enrichment,", ")[[1]]

##CIS 0hrs
CIS_0h_Med_ROS<-CIS_0h_cd[CIS_ROS_Med_genes,]

pcs <- MUDAN::getPcs(CIS_0h_Med_ROS,
                     nGenes=length(CIS_0h_Med_ROS), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X, Y=CIS_0h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df1, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 0h_ROS")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X, Y=CIS_0h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df1, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 0h_ROS")

##CIS 12hrs
CIS_12h_Med_ROS<-CIS_12h_cd[CIS_ROS_Med_genes,]

pcs <- MUDAN::getPcs(CIS_12h_Med_ROS,
                     nGenes=length(CIS_12h_Med_ROS), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X, Y=CIS_12h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df2, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 12h_ROS")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X, Y=CIS_12h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df2, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 12h_ROS")

##CIS 24hrs
CIS_24h_Med_ROS<-CIS_24h_cd[CIS_ROS_Med_genes,]

pcs <- MUDAN::getPcs(CIS_24h_Med_ROS,
                     nGenes=length(CIS_24h_Med_ROS), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X, Y=CIS_24h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df3, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 24h_ROS")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X, Y=CIS_24h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df3, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 24h_ROS")

##CIS 48hrs
CIS_48h_Med_ROS<-CIS_48h_cd[CIS_ROS_Med_genes,]

pcs <- MUDAN::getPcs(CIS_48h_Med_ROS,
                     nGenes=length(CIS_48h_Med_ROS), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X, Y=CIS_48h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df4, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 48h_ROS")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X, Y=CIS_48h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df4, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 48h_ROS")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TCA Cycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CIS_TCA_Med_genes<-strsplit(CIS_Medulla[grep("Citrate cycle (TCA cycle)", CIS_Medulla$Description ),]$core_enrichment,", ")[[1]]
CIS_TCA_Med_genes<-strsplit(CIS_Medulla[37,]$core_enrichment,", ")[[1]]

##CIS 0hrs
CIS_0h_Med_TCA<-CIS_0h_cd[CIS_TCA_Med_genes,]
set.seed(100)
pcs <- MUDAN::getPcs(CIS_0h_Med_TCA,
                     nGenes=length(CIS_0h_Med_TCA), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X, Y=CIS_0h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df1, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 0h_TCA")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z1<-PC1
Z1_max<-max(PC1)
Z1_min<-min(PC1)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X, Y=CIS_0h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df1, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 0h_TCA")

##CIS 12hrs
CIS_12h_Med_TCA<-CIS_12h_cd[CIS_TCA_Med_genes,]
set.seed(200)
pcs <- MUDAN::getPcs(CIS_12h_Med_TCA,
                     nGenes=length(CIS_12h_Med_TCA), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X, Y=CIS_12h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df2, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 12h_TCA")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z2<-PC1
Z2_max<-max(PC1)
Z2_min<-min(PC1)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X, Y=CIS_12h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df2, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 12h_TCA")

##CIS 24hrs
CIS_24h_Med_TCA<-CIS_24h_cd[CIS_TCA_Med_genes,]
set.seed(300)
pcs <- MUDAN::getPcs(CIS_24h_Med_TCA,
                     nGenes=length(CIS_24h_Med_TCA), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X, Y=CIS_24h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df3, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 24h_TCA")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z3<-PC1
Z3_max<-max(PC1)
Z3_min<-min(PC1)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X, Y=CIS_24h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df3, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 24h_TCA")

##CIS 48hrs
CIS_48h_Med_TCA<-CIS_48h_cd[CIS_TCA_Med_genes,]
set.seed(500)
pcs <- MUDAN::getPcs(CIS_48h_Med_TCA,
                     nGenes=length(CIS_48h_Med_TCA), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
PC1<-(PC1-min(PC1))/(max(PC1)-min(PC1))
names(PC1)<-rownames(pcs)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X, Y=CIS_48h$pos[names(PC1),]$Y, PC1=PC1)
ggplot(df4, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 48h_TCA")

#zscore of PC1
PC1<-data.frame(pcs)[,1]
PC1<-scale(PC1)
names(PC1)<-rownames(pcs)

Z4<-PC1
Z4_max<-max(PC1)
Z4_min<-min(PC1)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X, Y=CIS_48h$pos[names(PC1),]$Y, Zscore=PC1)
ggplot(df4, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 48h_TCA")

#Plot with common scale:
Zmax<-quantile(c(Z1,Z2,Z3,Z4), 0.99)
Zmin<-min(Z1_min, Z2_min, Z3_min, Z4_min)
Zmin

color_scale <- scale_color_gradientn(colors=viridis(5), values = c(0, 0.8, 1),
                                     limits = c(Zmin, Zmax), oob = scales::squish)

ggplot(df1, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 0h_TCA") + color_scale
ggplot(df2, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 12h_TCA") + color_scale
ggplot(df3, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 24h_TCA") + color_scale
ggplot(df4, aes(x=X, y=Y, col=Zscore)) + geom_point(size=0.5) + ggtitle("CIS 48h_TCA") + color_scale

################################################################################
#### PLOTTING: ASSUMING PC1 IS THE SAME FOR ALL TIMEPOINTS #####################
##CIS 0hrs
CIS_0h_Med_OXPHOS<-CIS_0h_cd[CIS_OXPHOS_Med_genes,]

pcs <- MUDAN::getPcs(CIS_0h_Med_OXPHOS,
                     nGenes=length(CIS_0h_Med_OXPHOS), 
                     nPcs=30,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
names(PC1)<-rownames(pcs)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X,
                Y=CIS_0h$pos[names(PC1),]$Y,
                PC1=PC1)

##CIS 12hrs
CIS_12h_Med_OXPHOS<-CIS_12h_cd[CIS_OXPHOS_Med_genes,]
pcs <- MUDAN::getPcs(CIS_12h_Med_OXPHOS,
                     nGenes=length(CIS_12h_Med_OXPHOS), 
                     nPcs=30,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
names(PC1)<-rownames(pcs)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X,
                Y=CIS_12h$pos[names(PC1),]$Y,
                PC1=PC1)

ggplot(df2, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 12h")

##CIS 24hrs
CIS_24h_Med_OXPHOS<-CIS_24h_cd[CIS_OXPHOS_Med_genes,]

pcs <- MUDAN::getPcs(CIS_24h_Med_OXPHOS,
                     nGenes=length(CIS_24h_Med_OXPHOS), 
                     nPcs=30,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
names(PC1)<-rownames(pcs)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X,
                Y=CIS_24h$pos[names(PC1),]$Y,
                PC1=PC1)
ggplot(df3, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 24h")

##CIS_48hrs
CIS_48h_Med_OXPHOS<-CIS_48h_cd[CIS_OXPHOS_Med_genes,]

pcs <- MUDAN::getPcs(CIS_48h_Med_OXPHOS,
                     nGenes=length(CIS_48h_Med_OXPHOS), 
                     nPcs=30,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
names(PC1)<-rownames(pcs)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X,
                Y=CIS_48h$pos[names(PC1),]$Y,
                PC1=PC1)
##Plotting
min<-min(df1$PC1, df2$PC2, df3$PC3, df4$PC1)
max<-max(df1$PC1, df2$PC2, df3$PC3, df4$PC1)

df1$PC1<-(df1$PC1-min)/(max-min)
df2$PC1<-(df2$PC1-min)/(max-min)
df3$PC1<-(df3$PC1-min)/(max-min)
df4$PC1<-(df4$PC1-min)/(max-min)

overall_min <- min(df1$PC1, df2$PC2, df3$PC3, df4$PC1) #for plot scale

overall_max <- quantile(c(df1$PC1, df2$PC2, df3$PC3, df4$PC1), 0.99) #for plot scale

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

p1<-ggplot(df1, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5, alpha=1) + color_scale +  ggtitle('CIS 0h')
p2<-ggplot(df2, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5, alpha=1) + color_scale +  ggtitle('CIS 12h') 
p3<-ggplot(df3, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5, alpha=1) + color_scale +  ggtitle('CIS 24h')
p4<-ggplot(df4, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5, alpha=1) + color_scale +  ggtitle('CIS 48h')

p1
p2
p3
p4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### GLYCOLYSIS
CIS_Glycolysis_Med_genes<-strsplit(CIS_Medulla[grep("Glycolysis / Gluconeogenesis", CIS_Medulla$Description ),]$core_enrichment,", ")[[1]]
CIS_Glycolysis_Med_genes<-strsplit(CIS_Medulla[grep("Oxidative phosphorylation", CIS_Medulla$Description ),]$core_enrichment,", ")[[1]]



##CIS 0hrs
CIS_0h_Med_Glycolysis<-CIS_0h_cd[CIS_Glycolysis_Med_genes,]

pcs <- MUDAN::getPcs(CIS_0h_Med_Glycolysis,
                     nGenes=length(CIS_0h_Med_Glycolysis), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
names(PC1)<-rownames(pcs)

df1<-data.frame(X=CIS_0h$pos[names(PC1),]$X,
                Y=CIS_0h$pos[names(PC1),]$Y,
                PC1=PC1)
ggplot(df1, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 0h")

##CIS 12hrs
CIS_12h_Med_Glycolysis<-CIS_12h_cd[CIS_Glycolysis_Med_genes,]

pcs <- MUDAN::getPcs(CIS_12h_Med_Glycolysis,
                     nGenes=length(CIS_12h_Med_Glycolysis), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
names(PC1)<-rownames(pcs)

df2<-data.frame(X=CIS_12h$pos[names(PC1),]$X,
                Y=CIS_12h$pos[names(PC1),]$Y,
                PC1=PC1)
ggplot(df2, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 12h")

##CIS 24hrs
CIS_24h_Med_Glycolysis<-CIS_24h_cd[CIS_Glycolysis_Med_genes,]

pcs <- MUDAN::getPcs(CIS_24h_Med_Glycolysis,
                     nGenes=length(CIS_24h_Med_Glycolysis), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
names(PC1)<-rownames(pcs)

df3<-data.frame(X=CIS_24h$pos[names(PC1),]$X,
                Y=CIS_24h$pos[names(PC1),]$Y,
                PC1=PC1)
ggplot(df3, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 24h")


##CIS 48hrs
CIS_48h_Med_Glycolysis<-CIS_48h_cd[CIS_Glycolysis_Med_genes,]

pcs <- MUDAN::getPcs(CIS_48h_Med_Glycolysis,
                     nGenes=length(CIS_48h_Med_Glycolysis), 
                     nPcs=5,
                     verbose=FALSE) 

PC1<-data.frame(pcs)[,1]
names(PC1)<-rownames(pcs)

df4<-data.frame(X=CIS_48h$pos[names(PC1),]$X,
                Y=CIS_48h$pos[names(PC1),]$Y,
                PC1=PC1)
ggplot(df4, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5) + ggtitle("CIS 48h")

##Plotting
min<-min(df1$PC1, df2$PC2, df3$PC3, df4$PC1)
max<-max(df1$PC1, df2$PC2, df3$PC3, df4$PC1)

df1$PC1<-(df1$PC1-min)/(max-min)
df2$PC1<-(df2$PC1-min)/(max-min)
df3$PC1<-(df3$PC1-min)/(max-min)
df4$PC1<-(df4$PC1-min)/(max-min)

overall_min <- min(df1$PC1, df2$PC2, df3$PC3, df4$PC1) #for plot scale

overall_max <- quantile(c(df1$PC1, df2$PC2, df3$PC3, df4$PC1), 0.99) #for plot scale

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

p1<-ggplot(df1, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5, alpha=1) + color_scale +  ggtitle('CIS 0h_Glycolysis')
p2<-ggplot(df2, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5, alpha=1) + color_scale +  ggtitle('CIS 12h_Glycolysis') 
p3<-ggplot(df3, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5, alpha=1) + color_scale +  ggtitle('CIS 24h_Glycolysis')
p4<-ggplot(df4, aes(x=X, y=Y, col=PC1)) + geom_point(size=0.5, alpha=1) + color_scale +  ggtitle('CIS 48h_Glycolysis')

p1
p2
p3
p4
