#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(enrichplot) #For ClusterProfiler()

############ KEGG GSEA PATHWAYS PLOTTING #######################################
########### CIS CORTEX (OXPHOS, Thermogenesis, TCA Cycle) ######################
#Extracting gene set enrichment analysis results (refer to the code 
#Linear_Regression_CIS_AKI.R stored in the the folder 2.Linear_Regression_Modeling)
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
#Extracting gene set enrichment analysis results (refer to the code 
#Linear_Regression_CIS_AKI.R stored in the the folder 2.Linear_Regression_Modeling)
#Note: INTERFACE is same as Outer Medulla (as described in the main article section)
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
#Extracting gene set enrichment analysis results (refer to the code 
#Linear_Regression_CIS_AKI.R stored in the the folder 2.Linear_Regression_Modeling)
#Note: MEDULLA is same as Inner Medulla (as described in the main article section)
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

