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
CIS_Cortex<-read.xlsx("Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Cortex")
CIS_Cortex$Description<-gsub(" - Mus.*","",CIS_Cortex$Description)

CIS_Interface<-read.xlsx("Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Interface")
CIS_Interface$Description<-gsub(" - Mus.*","",CIS_Interface$Description)

CIS_Medulla<-read.xlsx("Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Medulla")
CIS_Medulla$Description<-gsub(" - Mus.*","",CIS_Medulla$Description)

CIS_OXPHOS_Cor_genes<-strsplit(CIS_Cortex[grep("Oxidative phosphorylation", CIS_Cortex$Description ),]$core_enrichment,", ")[[1]]
CIS_OXPHOS_Int_genes<-strsplit(CIS_Interface[grep("Oxidative phosphorylation", CIS_Interface$Description ),]$core_enrichment,", ")[[1]]
CIS_OXPHOS_Med_genes<-strsplit(CIS_Medulla[grep("Oxidative phosphorylation", CIS_Medulla$Description ),]$core_enrichment,", ")[[1]]

Reduce(intersect, list(CIS_OXPHOS_Cor_genes, CIS_OXPHOS_Int_genes, CIS_OXPHOS_Med_genes))
intersect(CIS_OXPHOS_Cor_genes, CIS_OXPHOS_Med_genes)
intersect(CIS_OXPHOS_Cor_genes, CIS_OXPHOS_Int_genes)
intersect(CIS_OXPHOS_Med_genes, CIS_OXPHOS_Int_genes)
 