#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

## load data
load("data/AKI_data.RData")
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
library(cowplot)

################################################################################
########### CIS CORTEX (Hepatitis/C Epstein-Barr) ###########
CIS_Cortex.KEGG<-readRDS("EnrichmentPlots/CIS_Cortex.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(CIS_Cortex.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(CIS_Cortex.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_CIS_Cortex_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(CIS_Cortex.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="CIS_Cortex_Viral_Infection")
dev.off()

########### CIS INTERFACE (Hepatitis/C Epstein-Barr) ###########
CIS_Interface.KEGG<-readRDS("EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(CIS_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(CIS_Interface.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_CIS_Interface_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="CIS_Interface_Viral_Infection")
dev.off()

########### CIS MEDULLA (Hepatitis/C Epstein-Barr) ###########
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(CIS_Medulla.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_CIS_Medulla_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(CIS_Medulla.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="CIS_Medulla_Viral_Infection")
dev.off()

########### AKI CORTEX (Hepatitis/C Epstein-Barr) ###########
AKI_Cortex.KEGG<-readRDS("EnrichmentPlots/AKI_Cortex.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(AKI_Cortex.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(AKI_Cortex.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_AKI_Cortex_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(AKI_Cortex.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="AKI_Cortex_Viral_Infection")
dev.off()

########### AKI INTERFACE (Hepatitis/C Epstein-Barr) ###########
AKI_Interface.KEGG<-readRDS("EnrichmentPlots/AKI_Interface.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(AKI_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(AKI_Interface.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_AKI_Interface_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(AKI_Interface.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="AKI_Interface_Viral_Infection")
dev.off()

########### AKI MEDULLA (Hepatitis/C Epstein-Barr) ###########
AKI_Medulla.KEGG<-readRDS("EnrichmentPlots/AKI_Medulla.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(AKI_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(AKI_Medulla.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_AKI_Medulla_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(AKI_Medulla.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="AKI_Medulla_Viral_Infection")
dev.off()

########### CIS vs AKI INTERFACE (COMPLEMENT) ###########
CIS_Interface.HALL<-readRDS("EnrichmentPlots/CIS_Interface.HALL_GSEA.rds")
pathway1<-"COMPLEMENT"
ID1=which(CIS_Interface.HALL$Description==paste0("HALLMARK_", pathway1) )
pdf("Figures/Figure4/pdfs/GSEA_CIS_Interface_Complement.pdf", height=4, width=8)
gseaplot(CIS_Interface.KEGG, geneSetID = c(ID1),by="runningScore", title="CIS_Interface_COMPLEMENT")
dev.off()

AKI_Interface.HALL<-readRDS("EnrichmentPlots/AKI_Interface.HALL_GSEA.rds")
pathway1<-"COMPLEMENT"
ID1=which(AKI_Interface.HALL$Description==paste0("HALLMARK_", pathway1) )
pdf("Figures/Figure4/pdfs/GSEA_AKI_Interface_Complement.pdf", height=4, width=8)
gseaplot(AKI_Interface.HALL, geneSetID = c(ID1), by="runningScore", title="AKI_Interface_COMPLEMENT")
dev.off()

##### MITOPHAGY
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Mitophagy - animal"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pdf("Figures/FigureY/pdfs/GSEA_CIS_InnerMedulla_Mitophagy.pdf", height=4, width=8)
gseaplot(CIS_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="CIS_Imedulla_Mitophagy")
dev.off()
