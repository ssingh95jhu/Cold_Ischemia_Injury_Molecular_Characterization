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
library(cowplot)

############ KEGG GSEA PATHWAYS PLOTTING #######################################
########### CIS CORTEX (OXPHOS, Thermogenesis, TCA Cycle) ######################
CIS_Cortex.KEGG<-readRDS("EnrichmentPlots/CIS_Cortex.KEGG_GSEA.rds")

pathway1<-"Oxidative phosphorylation"
ID1=which(CIS_Cortex.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Citrate cycle (TCA cycle)"
ID2=which(CIS_Cortex.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pathway3<-"Thermogenesis"
ID3=which(CIS_Cortex.KEGG$Description==paste0(pathway3, " - Mus musculus (house mouse)") )

pdf("Figures/Figure3/pdfs/GSEA_Cortex_OXPHOS_TCA_Themogenesis.pdf", height=4, width=8)
gseaplot2(CIS_Cortex.KEGG, geneSetID = c(ID1,ID2,ID3),  subplots=1:2, color= c("brown","orange","red" ), title="CIS Cortex") 
dev.off()

########### CIS INTERFACE (OXPHOS, Thermogenesis, TCA Cycle)
CIS_Interface.KEGG<-readRDS("EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")

pathway1<-"Oxidative phosphorylation"
ID1=which(CIS_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Citrate cycle (TCA cycle)"
ID2=which(CIS_Interface.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pathway3<-"Thermogenesis"
ID3=which(CIS_Interface.KEGG$Description==paste0(pathway3, " - Mus musculus (house mouse)") )

pdf("Figures/Figure3/pdfs/GSEA_Interface_OXPHOS_TCA_Themogenesis.pdf", height=4, width=8)
gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1,ID2,ID3),  subplots=1:2, color= c("brown","orange","red" ), title="CIS Interface")
dev.off()

########### CIS MEDULLA (OXPHOS, Thermogenesis, TCA Cycle)
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")

pathway1<-"Oxidative phosphorylation"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Citrate cycle (TCA cycle)"
ID2=which(CIS_Medulla.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pathway3<-"Thermogenesis"
ID3=which(CIS_Medulla.KEGG$Description==paste0(pathway3, " - Mus musculus (house mouse)") )

pdf("Figures/Figure3/pdfs/GSEA_Medulla_OXPHOS_TCA_Themogenesis.pdf", height=4, width=8)
gseaplot2(CIS_Medulla.KEGG, geneSetID = c(ID1,ID2,ID3),  subplots=1:2, color= c("brown","orange","red" ), title="CIS Medulla")
dev.off()

