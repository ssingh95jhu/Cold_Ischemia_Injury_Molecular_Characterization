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
library(cowplot)
library(dplyr)
library(viridis) #color blind friendly
library(RColorBrewer) #for RColorbrewer

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

##############################################################################
#Extracting the compartment specific spots for CIS and AKI:
cortex <- readRDS('AKI-CIS-irl-ctrl_Cortex_spots.rds')
interface <- readRDS('AKI-CIS-irl-ctrl_Interface_spots.rds')
medulla <- readRDS('AKI-CIS-irl-ctrl_Medulla_spots.rds')

##x*****************************************************************************
#DESeq2 Operation to find top markers for each kidney compartment **************
##x*****************************************************************************
wb1<-createWorkbook()

all_spots<-data.frame(readRDS('AKI-CIS-irl-ctrl_Harmonized_AllClusters.rds'))
all_spots$Tissue<-'Other'
colnames(all_spots)<-c('Cluster','Tissue')
head(all_spots)

all_spots[names(cortex),]$Tissue<-'Cortex'
all_spots[names(interface),]$Tissue<-'Interface'
all_spots[names(medulla),]$Tissue<-'Medulla'

#Extracting all the CTRL spots:
CTRL_spots<-all_spots[grep("^CTRL", rownames(all_spots)),]
CTRL_spots.cortex<-CTRL_spots[which(CTRL_spots$Tissue=='Cortex'),]
CTRL_spots.interface<-CTRL_spots[which(CTRL_spots$Tissue=='Interface'),]
CTRL_spots.medulla<-CTRL_spots[which(CTRL_spots$Tissue=='Medulla'),]

#Extracting all the IRL (AKI24) spots:
IRL_spots<-all_spots[grep("^IR", rownames(all_spots)),]
IRL_spots.cortex<-IRL_spots[which(IRL_spots$Tissue=='Cortex'),]
IRL_spots.interface<-IRL_spots[which(IRL_spots$Tissue=='Interface'),]
IRL_spots.medulla<-IRL_spots[which(IRL_spots$Tissue=='Medulla'),]

# #Gene expression matrix from the all the kidney datasets
# cd<-cbind(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),],
#           CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),],
#           CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),],
#           CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),],
#           AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),],
#           AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),],
#           AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),],
#           AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),],
#           AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),],
#           ctrl1$gexp[unique(rownames(CIS_0h$gexp)),],
#           ctrl2$gexp[unique(rownames(CIS_0h$gexp)),],
#           ctrl3$gexp[unique(rownames(CIS_0h$gexp)),],
#           ctrl4$gexp[unique(rownames(CIS_0h$gexp)),],
#           irl1$gexp[unique(rownames(CIS_0h$gexp)),],
#           irl2$gexp[unique(rownames(CIS_0h$gexp)),],
#           irl3$gexp[unique(rownames(CIS_0h$gexp)),],
#           irl4$gexp[unique(rownames(CIS_0h$gexp)),] )

#Gene expression matrix from the all the kidney datasets:
cd<-cbind(ctrl1$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl2$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl3$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl4$gexp[unique(rownames(CIS_0h$gexp)),],
          irl1$gexp[unique(rownames(CIS_0h$gexp)),],
          irl2$gexp[unique(rownames(CIS_0h$gexp)),],
          irl3$gexp[unique(rownames(CIS_0h$gexp)),],
          irl4$gexp[unique(rownames(CIS_0h$gexp)),] )

#Subsetting the count matrix into different compartments:
CTRL.cd.cortex<-cd[,rownames(CTRL_spots.cortex)]
CTRL.cd.interface<-cd[,rownames(CTRL_spots.interface)]
CTRL.cd.medulla<-cd[,rownames(CTRL_spots.medulla)]

IRL.cd.cortex<-cd[,rownames(IRL_spots.cortex)]
IRL.cd.interface<-cd[,rownames(IRL_spots.interface)]
IRL.cd.medulla<-cd[,rownames(IRL_spots.medulla)]

#TOP MARKER GENES FOR EACH COMPARTMENT #########################################
#Performing DESeq2 Analysis between the CORTEX of CTRL vs IRL
set.seed(645)

rand1<-sample(colnames(CTRL.cd.cortex),1200)
rand2<-sample(colnames(IRL.cd.cortex),1200)

CTRL.cortex<-cd[,rand1]
dim(CTRL.cortex)

IRL.cortex<-cd[,rand2]
dim(IRL.cortex)

##x*****************************************************************************
#Preparing gene matrix for perfroming DESeq2 (Cortex)
compare.cd<-as.matrix(cbind(CTRL.cortex, IRL.cortex))
compare.cd<-compare.cd+1

cold<-cbind(c(rep("CTRL_Cortex",ncol(CTRL.cortex)), 
              rep("IRL_Cortex", ncol(IRL.cortex)) )) 

rownames(cold)<-c(colnames(compare.cd))
colnames(cold)<-c("condition")

#To check if all the genes (rows) are consistent across the different timepoints
all(rownames(cold) %in% colnames(compare.cd))
all(rownames(cold) == colnames(compare.cd))

#Creating DESeq2 object and performing differential expression analysis.
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = compare.cd, colData = cold, design = ~ condition)

#Defining the reference level for DESeq2 operation
#dds$condition <- relevel(dds$condition, ref = "0h") <-performing DESeq2 between all the conditions

keep <- rowSums(counts(dds)) > ncol(compare.cd) #Only genes with >=10 counts are considered
#keep <- rowSums(counts(dds)) > 10000
dds <- dds[keep,]

dds <- DESeq(dds)

#Saving the normalized count matrix
norm_counts<-counts(dds, normalized = TRUE)

#r<-results(dds)
res<-as.matrix(results(dds, contrast=c("condition", 'IRL_Cortex','CTRL_Cortex')) )

res1<-res[which(res[,6]<0.05),]
res1<-res1[rev(order(res1[,2])),]
res1<-cbind(rownames(res1),res1)
colnames(res1)[1]<-c('Genes')


#ovdisp_genes<-readRDS("OverDispersed_Features_for_PC.rds")
#cortex_markers<-intersect(rownames(res1), ovdisp_genes)
cortex_degs<-res1

# res2<-res1[cortex_markers,]
# cortex_markers<-res2[rev(order(res2[,3])),]
# head(res2)

addWorksheet(wb1, "CTRL_IRL_Cor")
writeData(wb1, "CTRL_IRL_Cor", cortex_degs[,c(1,3,7)])

##INTERFACE ********************************************************************
#Performing DESeq2 Analysis between the INTERFACE of CTRL vs IRL 
set.seed(805)

rand1<-sample(colnames(CTRL.cd.interface),1200)
rand2<-sample(colnames(IRL.cd.interface),1200)

CTRL.interface<-cd[,rand1]
dim(CTRL.interface)

IRL.interface<-cd[,rand2]
dim(IRL.interface)

##x*****************************************************************************
#Preparing gene matrix for perfroming DESeq2 (Interface)
compare.cd<-as.matrix(cbind(CTRL.interface, IRL.interface))
compare.cd<-compare.cd+1

cold<-cbind(c(rep("CTRL_Interface",ncol(CTRL.interface)), 
              rep("IRL_Interface", ncol(IRL.interface)) )) 

rownames(cold)<-c(colnames(compare.cd))
colnames(cold)<-c("condition")

#To check if all the genes (rows) are consistent across the different timepoints
all(rownames(cold) %in% colnames(compare.cd))
all(rownames(cold) == colnames(compare.cd))

#Creating DESeq2 object and performing differential expression analysis.
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = compare.cd, colData = cold, design = ~ condition)

#Defining the reference level for DESeq2 operation
#dds$condition <- relevel(dds$condition, ref = "0h") <-performing DESeq2 between all the conditions

keep <- rowSums(counts(dds)) > ncol(compare.cd) #Only genes with >=10 counts are considered
#keep <- rowSums(counts(dds)) > 10000
dds <- dds[keep,]

dds <- DESeq(dds)

#Saving the normalized count matrix
norm_counts<-counts(dds, normalized = TRUE)

#r<-results(dds)
res<-as.matrix(results(dds, contrast=c("condition", 'IRL_Interface','CTRL_Interface')) )

res1<-res[which(res[,6]<0.05),]
res1<-res1[rev(order(res1[,2])),]
res1<-cbind(rownames(res1),res1)
colnames(res1)[1]<-c('Genes')

interface_degs<-res1

# res2<-res1[cortex_markers,]
# cortex_markers<-res2[rev(order(res2[,3])),]
# head(res2)

addWorksheet(wb1, "CTRL_IRL_Int")
writeData(wb1, "CTRL_IRL_Int", interface_degs[,c(1,3,7)])


##MEDULLA **********************************************************************
#Performing DESeq2 Analysis between the MEDULLA of CTRL vs IRL
set.seed(355)

rand1<-sample(colnames(CTRL.cd.medulla),1200)
rand2<-sample(colnames(IRL.cd.medulla),1200)

CTRL.medulla<-cd[,rand1]
dim(CTRL.medulla)

IRL.medulla<-cd[,rand2]
dim(IRL.medulla)

##x*****************************************************************************
#Preparing gene matrix for perfroming DESeq2 (Interface)
compare.cd<-as.matrix(cbind(CTRL.medulla, IRL.medulla))
compare.cd<-compare.cd+1

cold<-cbind(c(rep("CTRL_Medulla",ncol(CTRL.medulla)), 
              rep("IRL_Medulla", ncol(IRL.medulla)) )) 

rownames(cold)<-c(colnames(compare.cd))
colnames(cold)<-c("condition")

#To check if all the genes (rows) are consistent across the different timepoints
all(rownames(cold) %in% colnames(compare.cd))
all(rownames(cold) == colnames(compare.cd))

#Creating DESeq2 object and performing differential expression analysis.
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = compare.cd, colData = cold, design = ~ condition)

#Defining the reference level for DESeq2 operation
#dds$condition <- relevel(dds$condition, ref = "0h") <-performing DESeq2 between all the conditions

keep <- rowSums(counts(dds)) > ncol(compare.cd) #Only genes with >=10 counts are considered
#keep <- rowSums(counts(dds)) > 10000
dds <- dds[keep,]

dds <- DESeq(dds)

#Saving the normalized count matrix
norm_counts<-counts(dds, normalized = TRUE)

#r<-results(dds)
res<-as.matrix(results(dds, contrast=c("condition", 'IRL_Medulla','CTRL_Medulla')) )

res1<-res[which(res[,6]<0.05),]
res1<-res1[rev(order(res1[,2])),]
res1<-cbind(rownames(res1),res1)
colnames(res1)[1]<-c('Genes')

medulla_degs<-res1
medulla_degs[,c(1,3,7)]

addWorksheet(wb1, "CTRL_IRL_Med")
writeData(wb1, "CTRL_IRL_Med", medulla_degs[,c(1,3,7)])

saveWorkbook(wb1, "Tables/Fig_4S_CTRL_vs_IRL_Compartmental_DEGs.xlsx", overwrite = TRUE)


