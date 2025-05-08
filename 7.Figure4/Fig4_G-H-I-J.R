#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(clusterProfiler)
library(enrichplot) #For ClusterProfiler()

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

########### CIS INTERFACE (Hepatitis/C Epstein-Barr) ###########################
#Note: INTERFACE is same as Outer Medulla (as described in the main article section)
CIS_Interface.KEGG<-readRDS("EnrichmentPlots/CIS_Interface.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(CIS_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(CIS_Interface.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_CIS_Interface_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(CIS_Interface.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="CIS_Interface_Viral_Infection")
dev.off()

########### CIS MEDULLA (Hepatitis/C Epstein-Barr) #############################
#Note: MEDULLA is same as Inner Medulla (as described in the main article section)
CIS_Medulla.KEGG<-readRDS("EnrichmentPlots/CIS_Medulla.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(CIS_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(CIS_Medulla.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_CIS_Medulla_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(CIS_Medulla.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="CIS_Medulla_Viral_Infection")
dev.off()

########### AKI CORTEX (Hepatitis/C Epstein-Barr) ##############################
AKI_Cortex.KEGG<-readRDS("EnrichmentPlots/AKI_Cortex.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(AKI_Cortex.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(AKI_Cortex.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_AKI_Cortex_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(AKI_Cortex.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="AKI_Cortex_Viral_Infection")
dev.off()

########### AKI INTERFACE (Hepatitis/C Epstein-Barr) ###########################
#Note: INTERFACE is same as Outer Medulla (as described in the main article section)
AKI_Interface.KEGG<-readRDS("EnrichmentPlots/AKI_Interface.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(AKI_Interface.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(AKI_Interface.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_AKI_Interface_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(AKI_Interface.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="AKI_Interface_Viral_Infection")
dev.off()

########### AKI MEDULLA (Hepatitis/C Epstein-Barr) #############################
#Note: MEDULLA is same as Inner Medulla (as described in the main article section)
AKI_Medulla.KEGG<-readRDS("EnrichmentPlots/AKI_Medulla.KEGG_GSEA.rds")
pathway1<-"Epstein-Barr virus infection"
ID1=which(AKI_Medulla.KEGG$Description==paste0(pathway1, " - Mus musculus (house mouse)") )

pathway2<-"Hepatitis C"
ID2=which(AKI_Medulla.KEGG$Description==paste0(pathway2, " - Mus musculus (house mouse)") )

pdf("Figures/Figure4/pdfs/GSEA_AKI_Medulla_Viral_HepatitisC_EpsteinBarr.pdf", height=4, width=8)
gseaplot2(AKI_Medulla.KEGG, geneSetID = c(ID1,ID2),  subplots=1:2, color= c("red", "black"), title="AKI_Medulla_Viral_Infection")
dev.off()

########### CIS vs AKI INTERFACE (COMPLEMENT) ##################################
#Note: INTERFACE is same as Outer Medulla (as described in the main article section)
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
