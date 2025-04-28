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

########### CIS vs AKI  (Glycolysis in Medulla) ##################################
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Glycolysis / Gluconeogenesis"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)" ))

pdf("Figures/Figure4S/pdfs/GSEA_CIS_Medulla_Glycolysis.pdf", height=4, width=8)
gseaplot(CIS_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="CIS_Medulla_Glycolysis")
dev.off()

AKI_Medulla.KEGG<-readRDS("EnrichmentPlots/AKI_Medulla.KEGG_GSEA.rds")
pathway1<-"Glycolysis / Gluconeogenesis"
ID1=which(AKI_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)" ))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_Medulla_Glycolysis.pdf", height=4, width=8)
gseaplot(AKI_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="AKI_Medulla_Glycolysis")
dev.off()

########### CIS vs AKI  (Oxocarboxylic in Medulla) ##################################
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"2-Oxocarboxylic acid metabolism"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)" ))

pdf("Figures/Figure4S/pdfs/GSEA_CIS_Medulla_OxocarboxylicAcid.pdf", height=4, width=8)
gseaplot(CIS_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="CIS_Medulla_2OxocarboxylicAcid")
dev.off()

AKI_Medulla.KEGG<-readRDS("EnrichmentPlots/AKI_Medulla.KEGG_GSEA.rds")
pathway1<-"2-Oxocarboxylic acid metabolism"
ID1=which(AKI_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)" ))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_Medulla_OxocarboxylicAcid.pdf", height=4, width=8)
gseaplot(AKI_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="AKI_Medulla_2OxocarboxylicAcid")
dev.off()

########### CIS vs AKI  (TCA Cycle in Medulla) ##################################
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Citrate cycle (TCA cycle)"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)" ))

pdf("Figures/Figure4S/pdfs/GSEA_CIS_Medulla_TCAcycle.pdf", height=4, width=8)
gseaplot(CIS_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="CIS_Medulla_TCAcycle")
dev.off()

AKI_Medulla.KEGG<-readRDS("EnrichmentPlots/AKI_Medulla.KEGG_GSEA.rds")
pathway1<-"Citrate cycle (TCA cycle)"
ID1=which(AKI_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)" ))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_Medulla_TCAcycle.pdf", height=4, width=8)
gseaplot(AKI_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="AKI_Medulla_TCAcycle")
dev.off()

########### CIS vs AKI  (OXPHOS in Medulla) ##################################
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Oxidative phosphorylation"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)" ))

pdf("Figures/Figure4S/pdfs/GSEA_CIS_Medulla_OXPHOS.pdf", height=4, width=8)
gseaplot(CIS_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="CIS_Medulla_OXPHOS")
dev.off()

AKI_Medulla.KEGG<-readRDS("EnrichmentPlots/AKI_Medulla.KEGG_GSEA.rds")
pathway1<-"Oxidative phosphorylation"
ID1=which(AKI_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)" ))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_Medulla_OXPHOS.pdf", height=4, width=8)
gseaplot(AKI_Medulla.KEGG, geneSetID = c(ID1),by="runningScore", title="AKI_Medulla_OXPHOS")
dev.off()