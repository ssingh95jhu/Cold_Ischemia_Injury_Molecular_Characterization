#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

## load data
load("data/CIS_data.RData")
load("data/AKI_data.RData")
load("data/Rabb_ctrl.RData")
load("data/Rabb_irl24h.RData")

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(cowplot)
library(viridis) #color blind friendly
library(RColorBrewer) #for RColorbrewer

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

#################### NORMALIZE IRL DATASET ####################################

## limma performs better with gaussian distributed
irl1$mat <- MERINGUE::normalizeCounts(irl1$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
irl2$mat <- MERINGUE::normalizeCounts(irl2$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
irl3$mat <- MERINGUE::normalizeCounts(irl3$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
irl4$mat <- MERINGUE::normalizeCounts(irl4$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)

## easier for plotting
irl1$mat_notlog <- MERINGUE::normalizeCounts(irl1$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl2$mat_notlog <- MERINGUE::normalizeCounts(irl2$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl3$mat_notlog <- MERINGUE::normalizeCounts(irl3$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl4$mat_notlog <- MERINGUE::normalizeCounts(irl4$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

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

##############################################################################
#Extracting the compartment specific spots for CIS and AKI:

cortex <- readRDS('AKI-CIS-irl-ctrl_Cortex_spots.rds')
interface <- readRDS('AKI-CIS-irl-ctrl_Interface_spots.rds')
medulla <- readRDS('AKI-CIS-irl-ctrl_Medulla_spots.rds')
#all_spots<-readRDS('AKI-CIS-irl-ctrl_All_spots.rds')

#1.All_Cortex_Spots (CTRL)
CTRL1_cortex<-names(cortex[grep("^CTRL1", names(cortex))])
CTRL2_cortex<-names(cortex[grep("^CTRL2", names(cortex))])
CTRL3_cortex<-names(cortex[grep("^CTRL3", names(cortex))])
CTRL4_cortex<-names(cortex[grep("^CTRL4", names(cortex))])

#2.All_Interface_Spots (CTRL)
CTRL1_interface<-names(interface[grep("^CTRL1", names(interface))])
CTRL2_interface<-names(interface[grep("^CTRL2", names(interface))])
CTRL3_interface<-names(interface[grep("^CTRL3", names(interface))])
CTRL4_interface<-names(interface[grep("^CTRL4", names(interface))])

#3.All_Medulla_Spots (CTRL)
CTRL1_medulla<-names(medulla[grep("^CTRL1", names(medulla))])
CTRL2_medulla<-names(medulla[grep("^CTRL2", names(medulla))])
CTRL3_medulla<-names(medulla[grep("^CTRL3", names(medulla))])
CTRL4_medulla<-names(medulla[grep("^CTRL4", names(medulla))])

#1.All_Cortex_Spots (IRL)
IRL1_cortex<-names(cortex[grep("^IR1", names(cortex))])
IRL2_cortex<-names(cortex[grep("^IR2", names(cortex))])
IRL3_cortex<-names(cortex[grep("^IR3", names(cortex))])
IRL4_cortex<-names(cortex[grep("^IR4", names(cortex))])

#2.All_Interface_Spots (IRL)
IRL1_interface<-names(interface[grep("^IR1", names(interface))])
IRL2_interface<-names(interface[grep("^IR2", names(interface))])
IRL3_interface<-names(interface[grep("^IR3", names(interface))])
IRL4_interface<-names(interface[grep("^IR4", names(interface))])

#3.All_Medulla_Spots (IRL)
IRL1_medulla<-names(medulla[grep("^IR1", names(medulla))])
IRL2_medulla<-names(medulla[grep("^IR2", names(medulla))])
IRL3_medulla<-names(medulla[grep("^IR3", names(medulla))])
IRL4_medulla<-names(medulla[grep("^IR4", names(medulla))])


#1.All_Cortex_Spots (CIS)
CIS_cortex.0h<-names(cortex[grep("^CIS_0h", names(cortex))])
CIS_cortex.12h<-names(cortex[grep("^CIS_12h", names(cortex))])
CIS_cortex.24h<-names(cortex[grep("^CIS_24h", names(cortex))])
CIS_cortex.48h<-names(cortex[grep("^CIS_48h", names(cortex))])

#2.All_Interface_Spots (CIS)
CIS_interface.0h<-names(interface[grep("^CIS_0h", names(interface))])
CIS_interface.12h<-names(interface[grep("^CIS_12h", names(interface))])
CIS_interface.24h<-names(interface[grep("^CIS_24h", names(interface))])
CIS_interface.48h<-names(interface[grep("^CIS_48h", names(interface))])

#3.All_Medulla_Spots (CIS)
CIS_medulla.0h<-names(medulla[grep("^CIS_0h", names(medulla))])
CIS_medulla.12h<-names(medulla[grep("^CIS_12h", names(medulla))])
CIS_medulla.24h<-names(medulla[grep("^CIS_24h", names(medulla))])
CIS_medulla.48h<-names(medulla[grep("^CIS_48h", names(medulla))])


#1.All_Cortex_Spots (AKI)
AKI_cortex.sham<-names(cortex[grep("^AKI_sham", names(cortex))])
AKI_cortex.4h<-names(cortex[grep("^AKI_4h", names(cortex))])
AKI_cortex.12h<-names(cortex[grep("^AKI_12h", names(cortex))])
AKI_cortex.2d<-names(cortex[grep("^AKI_2d", names(cortex))])
AKI_cortex.6w<-names(cortex[grep("^AKI_6w", names(cortex))])

#2.All_Interface_Spots (AKI)
AKI_interface.sham<-names(interface[grep("^AKI_sham", names(interface))])
AKI_interface.4h<-names(interface[grep("^AKI_4h", names(interface))])
AKI_interface.12h<-names(interface[grep("^AKI_12h", names(interface))])
AKI_interface.2d<-names(interface[grep("^AKI_2d", names(interface))])
AKI_interface.6w<-names(interface[grep("^AKI_6w", names(interface))])

#3.All_Medulla_Spots (AKI)
AKI_medulla.sham<-names(medulla[grep("^AKI_sham", names(medulla))])
AKI_medulla.4h<-names(medulla[grep("^AKI_4h", names(medulla))])
AKI_medulla.12h<-names(medulla[grep("^AKI_12h", names(medulla))])
AKI_medulla.2d<-names(medulla[grep("^AKI_2d", names(medulla))])
AKI_medulla.6w<-names(medulla[grep("^AKI_6w", names(medulla))])

####################### PLOTTING AND VISUALIZATION #############################




##################### PLOTTING (CTRL vs CIS vs AKI vs IRL) ############################
                     