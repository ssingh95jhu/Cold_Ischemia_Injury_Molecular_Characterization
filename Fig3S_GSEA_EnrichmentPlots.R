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

################################################################################
########### CIS CORTEX (Glycolysis, Fatty Acid, Tryptophan, Leucine) ###########
CIS_Cortex.KEGG<-readRDS("EnrichmentPlots/CIS_Cortex.KEGG_GSEA.rds")
pathway1<-"Fatty acid degradation"
ID1=which(CIS_Cortex.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Tryptophan metabolism"
ID2=which(CIS_Cortex.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pathway3<-"Valine, leucine and isoleucine degradation"
ID3=which(CIS_Cortex.KEGG$Description==paste0(pathway3, " - Mus musculus (house mouse)") )

pdf("Figures/Figure3S/pdfs/GSEA_Cortex_Glycolysis_FattyAcid_Tryptophan_Valine.pdf", height=4, width=8)
gseaplot2(CIS_Cortex.KEGG, geneSetID = c(ID1,ID2,ID3),  subplots=1:2, color= c("cyan","skyblue", "navy" ), title="CIS_Cortex_Metabolic_Substrate")
dev.off()


CIS_Interface.KEGG<-readRDS("EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")
pathway1<-"Fatty acid degradation"
ID1=which(CIS_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Tryptophan metabolism"
ID2=which(CIS_Interface.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pathway3<-"Valine, leucine and isoleucine degradation"
ID3=which(CIS_Interface.KEGG$Description==paste0(pathway3, " - Mus musculus (house mouse)") )

pdf("Figures/Figure3S/pdfs/GSEA_Interface_Glycolysis_FattyAcid_Tryptophan_Valine.pdf", height=4, width=8)
gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1,ID2,ID3),  subplots=1:2, color= c("cyan","skyblue", "navy"), title="CIS_Interface_Metabolic_Substrate")
dev.off()


CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Fatty acid degradation"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Tryptophan metabolism"
ID2=which(CIS_Medulla.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pathway3<-"Valine, leucine and isoleucine degradation"
ID3=which(CIS_Medulla.KEGG$Description==paste0(pathway3, " - Mus musculus (house mouse)") )

pathway4<-"Glycolysis / Gluconeogenesis"
ID4=which(CIS_Medulla.KEGG$Description==paste0(pathway4, " - Mus musculus (house mouse)") )

pdf("Figures/Figure3S/pdfs/GSEA_Medulla_Glycolysis_FattyAcid_Tryptophan_Valine.pdf", height=4, width=8)
gseaplot2(CIS_Medulla.KEGG, geneSetID = c(ID1,ID2,ID3, ID4), subplots=1:2, color= c("cyan","royalblue" ,"skyblue", "navy"), title="CIS_Medulla_Metabolic_Substrate")
dev.off()

########### CIS CORTEX (ROS) ###################################################
CIS_Cortex.KEGG<-readRDS("EnrichmentPlots/CIS_Cortex.KEGG_GSEA.rds")
pathway1<-"Chemical carcinogenesis - reactive oxygen species"
ID1=which(CIS_Cortex.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )
pdf("Figures/Figure3S/pdfs/GSEA_Cortex_ROS.pdf", height=4, width=8)
gseaplot(CIS_Cortex.KEGG, geneSetID = c(ID1),  by="runningScore", title="CIS_Cortex_ROS") 
#gseaplot2(CIS_Cortex.KEGG, geneSetID = c(ID1), subplots=1:2,  title="CIS Cortex_ROS")
dev.off()

CIS_Interface.KEGG<-readRDS("EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")
pathway1<-"Chemical carcinogenesis - reactive oxygen species"
ID1=which(CIS_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )
pdf("Figures/Figure3S/pdfs/GSEA_Interface_ROS.pdf", height=4, width=8)
gseaplot(CIS_Interface.KEGG, geneSetID = c(ID1),  by="runningScore", title="CIS_Interface_ROS") 
#gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1), subplots=1:2, title="CIS_Interface_ROS")
dev.off()

CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Chemical carcinogenesis - reactive oxygen species"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )
pdf("Figures/Figure3S/pdfs/GSEA_Medulla_ROS.pdf", height=4, width=8)
gseaplot(CIS_Medulla.KEGG, geneSetID = c(ID1),  by="runningScore", title="CIS_Medulla_ROS")  
#gseaplot2(CIS_Medulla.KEGG, geneSetID = c(ID1),  subplots=1:2, color= c("black"), title="CIS_Medulla_ROS") 
dev.off()

########### KEGG GSEA PATHWAYS PLOTTING ########################################
########### CIS CORTEX (Antigen Presentation) ##################################
CIS_Cortex.KEGG<-readRDS("EnrichmentPlots/CIS_Cortex.KEGG_GSEA.rds")
pathway1<-"Antigen processing and presentation"
ID1=which(CIS_Cortex.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )
pdf("Figures/Figure3S/pdfs/GSEA_Cortex_Antigen_Processing_and_Presentation.pdf.pdf", height=4, width=8)
gseaplot(CIS_Cortex.KEGG, geneSetID = c(ID1),  by="runningScore", title="CIS Cortex_Antigen_Presentation")  
#gseaplot2(CIS_Cortex.KEGG, geneSetID = c(ID1),  subplots=1:2, color= c("black"), title="CIS_Cortex_Antigen_Presentation")  
dev.off()

CIS_Interface.KEGG<-readRDS("EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")
pathway1<-"Antigen processing and presentation"
ID1=which(CIS_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )
pdf("Figures/Figure3S/pdfs/GSEA_Interface_Antigen_Processing_and_Presentation.pdf.pdf", height=4, width=8)
gseaplot(CIS_Interface.KEGG, geneSetID = c(ID1),  by="runningScore", title="CIS_Interface_Antigen_Presentation")  
#gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1),  subplots=1:2, color= c("black"), title="CIS_Interface_Antigen_Presentation") 
dev.off()

CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Antigen processing and presentation"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )
pdf("Figures/Figure3S/pdfs/GSEA_Medulla_Antigen_Processing_and_Presentation.pdf", height=4, width=8)
gseaplot(CIS_Medulla.KEGG, geneSetID = c(ID1),  by="runningScore", title="CIS_Medulla_Antigen_Presentation")  
#gseaplot2(CIS_Medulla.KEGG, geneSetID = c(ID1),  subplots=1:2, color= c("black"), title="CIS_Medulla_Antigen_Presentation") 
dev.off()