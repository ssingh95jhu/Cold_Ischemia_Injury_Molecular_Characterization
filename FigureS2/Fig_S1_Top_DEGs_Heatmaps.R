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

#Gene expression matrix from the all the kidney datasets
cd<-cbind(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),],
          CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),],
          CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),],
          CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl1$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl2$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl3$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl4$gexp[unique(rownames(CIS_0h$gexp)),],
          irl1$gexp[unique(rownames(CIS_0h$gexp)),],
          irl2$gexp[unique(rownames(CIS_0h$gexp)),],
          irl3$gexp[unique(rownames(CIS_0h$gexp)),],
          irl4$gexp[unique(rownames(CIS_0h$gexp)),] )

#Subsetting the count matrix into different compartments:
cd.cortex<-cd[,all_spots$Tissue=='Cortex']
cd.interface<-cd[,all_spots$Tissue=='Interface']
cd.medulla<-cd[,all_spots$Tissue=='Medulla']
cd.other<-cd[,all_spots$Tissue=='Other']

#TOP MARKER GENES FOR EACH COMPARTMENT #########################################
#Performing DESeq2 Analysis between the CORTEX vs DIFFERENT COMPARTMENTS
set.seed(245)
rand1<-sample(rownames(all_spots[all_spots$Tissue=='Cortex',]),2000)
rand2<-sample(rownames(all_spots[all_spots$Tissue!='Cortex',]),2000)

cd.interest<-cd[,rand1]
dim(cd.interest)
cd.noninterest<-cd[,rand2]

##x*****************************************************************************
#Preparing gene matrix for perfroming DESeq2 (Cortex)
compare.cd<-as.matrix(cbind(cd.interest, cd.noninterest))
compare.cd<-compare.cd+1

cold<-cbind(c(rep("Cortex",ncol(cd.interest)), 
              rep("NotCortex", ncol(cd.noninterest)) )) 

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
res<-as.matrix(results(dds, contrast=c("condition", 'Cortex','NotCortex')) )
res1<-res[which(res[,6]<0.05),]
res1<-res1[rev(order(res1[,2])),]
res1<-cbind(rownames(res1),res1)
colnames(res1)[1]<-c('Genes')

ovdisp_genes<-readRDS("OverDispersed_Features_for_PC.rds")
cortex_markers<-intersect(rownames(res1), ovdisp_genes)

res2<-res1[cortex_markers,]
cortex_markers<-res2[rev(order(res2[,3])),]
head(res2)

addWorksheet(wb1, "Cortex")
writeData(wb1, "Cortex", cortex_markers)

##INTERFACE ********************************************************************
#Performing DESeq2 Analysis between the INTERFACE vs DIFFERENT COMPARTMENTS
set.seed(255)
rand1<-sample(rownames(all_spots[all_spots$Tissue=='Interface',]),2000)
rand2<-sample(rownames(all_spots[all_spots$Tissue!='Interface',]),2000)

cd.interest<-cd[,rand1]
dim(cd.interest)
cd.noninterest<-cd[,rand2]

##x*****************************************************************************
#Preparing gene matrix for perfroming DESeq2 (Interface)
compare.cd<-as.matrix(cbind(cd.interest, cd.noninterest))
compare.cd<-compare.cd+1

cold<-cbind(c(rep("Interface",ncol(cd.interest)), 
              rep("NotInterface", ncol(cd.noninterest)) )) 

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
res<-as.matrix(results(dds, contrast=c("condition", 'Interface','NotInterface')) )
res1<-res[which(res[,6]<0.05),]
res1<-res1[rev(order(res1[,2])),]
res1<-cbind(rownames(res1),res1)
colnames(res1)[1]<-c('Genes')

#ovdisp_genes<-readRDS("OverDispersed_Features_for_PC.rds")
interface_markers<-intersect(rownames(res1), ovdisp_genes)

res2<-res1[interface_markers,]
interface_markers<-res2[rev(order(res2[,3])),]
head(interface_markers)

addWorksheet(wb1, "Interface")
writeData(wb1, "Interface", interface_markers)

##MEDULLA **********************************************************************
#Performing DESeq2 Analysis between the MEDULLA vs DIFFERENT COMPARTMENTS
set.seed(265)
rand1<-sample(rownames(all_spots[all_spots$Tissue=='Medulla',]),2000)
rand2<-sample(rownames(all_spots[all_spots$Tissue!='Medulla',]),2000)

cd.interest<-cd[,rand1]
dim(cd.interest)
cd.noninterest<-cd[,rand2]

##x*****************************************************************************
#Preparing gene matrix for perfroming DESeq2 (Medulla)
compare.cd<-as.matrix(cbind(cd.interest, cd.noninterest))
compare.cd<-compare.cd+1

cold<-cbind(c(rep("Medulla",ncol(cd.interest)), 
              rep("NotMedulla", ncol(cd.noninterest)) )) 

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
res<-as.matrix(results(dds, contrast=c("condition", 'Medulla','NotMedulla')) )
res1<-res[which(res[,6]<0.05),]
res1<-res1[rev(order(res1[,2])),]
res1<-cbind(rownames(res1),res1)
colnames(res1)[1]<-c('Genes')

#ovdisp_genes<-readRDS("OverDispersed_Features_for_PC.rds")
medulla_markers<-intersect(rownames(res1), ovdisp_genes)

res2<-res1[medulla_markers,]
medulla_markers<-res2[rev(order(res2[,3])),]
head(res2)

addWorksheet(wb1, "Medulla")
writeData(wb1, "Medulla", medulla_markers)

##OTHER **********************************************************************
#Performing DESeq2 Analysis between the OTHER vs DIFFERENT COMPARTMENTS
set.seed(275)
rand1<-sample(rownames(all_spots[all_spots$Tissue=='Other',]),2000)
rand2<-sample(rownames(all_spots[all_spots$Tissue!='Other',]),2000)

cd.interest<-cd[,rand1]
dim(cd.interest)
cd.noninterest<-cd[,rand2]

##x*****************************************************************************
#Preparing gene matrix for perfroming DESeq2 (Other)
compare.cd<-as.matrix(cbind(cd.interest, cd.noninterest))
compare.cd<-compare.cd+1

cold<-cbind(c(rep("Other",ncol(cd.interest)), 
              rep("NotOther", ncol(cd.noninterest)) )) 

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
res<-as.matrix(results(dds, contrast=c("condition", 'Other','NotOther')) )
res1<-res[which(res[,6]<0.05),]
res1<-res1[rev(order(res1[,2])),]
res1<-cbind(rownames(res1),res1)
colnames(res1)[1]<-c('Genes')

#ovdisp_genes<-readRDS("OverDispersed_Features_for_PC.rds")
other_markers<-intersect(rownames(res1), ovdisp_genes)

res2<-res1[other_markers,]
other_markers<-res2[rev(order(res2[,3])),]
head(res2)

addWorksheet(wb1, "Other")
writeData(wb1, "Other", other_markers)

saveWorkbook(wb1, "Tables/TopMarkers_Kidney_Compartments.xlsx")


############################ normalize #########################################
## limma performs better with gaussian distributed
## easier for plotting
CIS_0h$mat_notlog <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_12h$mat_notlog <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_24h$mat_notlog <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_48h$mat_notlog <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_cd<-cbind(CIS_0h$mat_notlog, CIS_12h$mat_notlog, CIS_24h$mat_notlog, CIS_48h$mat_notlog)

AKI_sham$mat_notlog <- MERINGUE::normalizeCounts(AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_4h$mat_notlog <- MERINGUE::normalizeCounts(AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_12h$mat_notlog <- MERINGUE::normalizeCounts(AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_2d$mat_notlog <- MERINGUE::normalizeCounts(AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_6w$mat_notlog <- MERINGUE::normalizeCounts(AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_cd<-cbind(AKI_sham$mat_notlog, AKI_4h$mat_notlog, AKI_12h$mat_notlog, AKI_2d$mat_notlog, AKI_6w$mat_notlog)

ctrl1$mat_notlog <- MERINGUE::normalizeCounts(ctrl1$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl2$mat_notlog <- MERINGUE::normalizeCounts(ctrl2$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl3$mat_notlog <- MERINGUE::normalizeCounts(ctrl3$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl4$mat_notlog <- MERINGUE::normalizeCounts(ctrl4$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl_cd<-cbind(ctrl1$mat_notlog, ctrl2$mat_notlog, ctrl3$mat_notlog, ctrl4$mat_notlog)

irl1$mat_notlog <- MERINGUE::normalizeCounts(irl1$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl2$mat_notlog <- MERINGUE::normalizeCounts(irl2$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl3$mat_notlog <- MERINGUE::normalizeCounts(irl3$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl4$mat_notlog <- MERINGUE::normalizeCounts(irl4$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl_cd<-cbind(irl1$mat_notlog, irl2$mat_notlog, irl3$mat_notlog, irl4$mat_notlog)


# ############################## CIS DATASET SPOTS ###############################
# #1.All_Cortex_Spots (CIS)
# CIS_cortex.0h<-names(cortex[grep("^CIS_0h", names(cortex))])
# CIS_cortex.12h<-names(cortex[grep("^CIS_12h", names(cortex))])
# CIS_cortex.24h<-names(cortex[grep("^CIS_24h", names(cortex))])
# CIS_cortex.48h<-names(cortex[grep("^CIS_48h", names(cortex))])
# 
# #2.All_Interface_Spots (CIS)
# CIS_interface.0h<-names(interface[grep("^CIS_0h", names(interface))])
# CIS_interface.12h<-names(interface[grep("^CIS_12h", names(interface))])
# CIS_interface.24h<-names(interface[grep("^CIS_24h", names(interface))])
# CIS_interface.48h<-names(interface[grep("^CIS_48h", names(interface))])
# 
# #3.All_Medulla_Spots (CIS)
# CIS_medulla.0h<-names(medulla[grep("^CIS_0h", names(medulla))])
# CIS_medulla.12h<-names(medulla[grep("^CIS_12h", names(medulla))])
# CIS_medulla.24h<-names(medulla[grep("^CIS_24h", names(medulla))])
# CIS_medulla.48h<-names(medulla[grep("^CIS_48h", names(medulla))])
# 
# ############################## AKI DATASET SPOTS ###############################
# #1.All_Cortex_Spots (AKI)
# AKI_cortex.sham<-names(cortex[grep("^AKI_sham", names(cortex))])
# AKI_cortex.4h<-names(cortex[grep("^AKI_4h", names(cortex))])
# AKI_cortex.12h<-names(cortex[grep("^AKI_12h", names(cortex))])
# AKI_cortex.2d<-names(cortex[grep("^AKI_2d", names(cortex))])
# 
# #2.All_Interface_Spots (AKI)
# AKI_interface.sham<-names(interface[grep("^AKI_sham", names(interface))])
# AKI_interface.4h<-names(interface[grep("^AKI_4h", names(interface))])
# AKI_interface.12h<-names(interface[grep("^AKI_12h", names(interface))])
# AKI_interface.2d<-names(interface[grep("^AKI_2d", names(interface))])
# 
# #3.All_Medulla_Spots (AKI)
# AKI_medulla.sham<-names(medulla[grep("^AKI_sham", names(medulla))])
# AKI_medulla.4h<-names(medulla[grep("^AKI_4h", names(medulla))])
# AKI_medulla.12h<-names(medulla[grep("^AKI_12h", names(medulla))])
# AKI_medulla.2d<-names(medulla[grep("^AKI_2d", names(medulla))])
# 
# 
# 
# CIS_cortex.0h.cd<-as.matrix(CIS_0h$mat_notlog[,CIS_cortex.0h])
# CIS_cortex.12h.cd<-as.matrix(CIS_12h$mat_notlog[,CIS_cortex.12h])
# CIS_interface.0h.cd<-as.matrix(CIS_0h$mat_notlog[,CIS_interface.0h])
# CIS_interface.12h.cd<-as.matrix(CIS_12h$mat_notlog[,CIS_interface.12h])
# CIS_medulla.0h.cd<-as.matrix(CIS_0h$mat_notlog[,CIS_medulla.0h])
# CIS_medulla.12h.cd<-as.matrix(CIS_12h$mat_notlog[,CIS_medulla.12h])
# 
# # Create a data frame
# g<-'Lcn2'
# data <- data.frame(
#   value = c(CIS_cortex.0h.cd[g,], CIS_cortex.12h.cd[g,],
#             CIS_interface.0h.cd[g,], CIS_interface.12h.cd[g,],
#             CIS_medulla.0h.cd[g,], CIS_medulla.12h.cd[g,]),
#   group = factor(c(rep("CIS_Cortex_0h", length(CIS_cortex.0h.cd[g,])),
#                  rep("CIS_Cortex_12h", length(CIS_cortex.12h.cd[g,])),
#                  rep("CIS_Interface_0h", length(CIS_interface.0h.cd[g,])),
#                  rep("CIS_Interface_12h", length(CIS_interface.12h.cd[g,])),
#                  rep("CIS_Medulla_0h", length(CIS_medulla.0h.cd[g,])),
#                  rep("CIS_Medulla_12h", length(CIS_medulla.12h.cd[g,]))
#                  ) ) )
# 
# # Load ggplot2
# library(ggplot2)
# 
# # Plot the violin plot
# ggplot(data, aes(x = group, y = value, fill = group)) +
#   geom_violin(trim = FALSE, alpha = 0.7, width=1) +  # Violin plot with transparency
#   geom_boxplot(width = 0.01, outlier.shape = NA,   # Boxplot with reduced width
#                alpha = 0.9, color = "black") +    # Boxplot with black outline
#   
#   theme_minimal() +          # Use a clean theme
#   labs(x = "Group", y = "Value") +
#   scale_fill_brewer(palette = "Set3") + # Optional color palette
#   theme(legend.position = "none") 
# 
# range(data$value)

################################################################################
################### TYPe 2 PLOT ################################################
#1.CIS Spots Count Matrix (CPM)
CIS_cortex.cd<-as.matrix(CIS_cd[,names(cortex[grep("^CIS", names(cortex))]) ])
CIS_interface.cd<-as.matrix(CIS_cd[,names(interface[grep("^CIS", names(interface))]) ])
CIS_medulla.cd<-as.matrix(CIS_cd[,names(medulla[grep("^CIS", names(medulla))]) ])

#2.AKI Spots Count Matrix (CPM)
AKI_cortex.cd<-as.matrix(AKI_cd[,names(cortex[grep("^AKI", names(cortex))]) ])
AKI_interface.cd<-as.matrix(AKI_cd[,names(interface[grep("^AKI", names(interface))]) ])
AKI_medulla.cd<-as.matrix(AKI_cd[,names(medulla[grep("^AKI", names(medulla))]) ])

#3.CTRL Spots Count Matrix (CPM)
CTRL_cortex.cd<-as.matrix(ctrl_cd[,names(cortex[grep("^CTRL", names(cortex))]) ])
CTRL_interface.cd<-as.matrix(ctrl_cd[,names(interface[grep("^CTRL", names(interface))]) ])
CTRL_medulla.cd<-as.matrix(ctrl_cd[,names(medulla[grep("^CTRL", names(medulla))]) ])

#4.IRL Spots Count Matrix (CPM)
IRL_cortex.cd<-as.matrix(irl_cd[,names(cortex[grep("^IR", names(cortex))]) ])
IRL_interface.cd<-as.matrix(irl_cd[,names(interface[grep("^IR", names(interface))]) ])
IRL_medulla.cd<-as.matrix(irl_cd[,names(medulla[grep("^IR", names(medulla))]) ])


#Make violin plot for the gene 'g'
g<-'Slc14a2'
data <- data.frame(
  value = c(log10(CIS_cortex.cd[g,]+1)/max(log10(CIS_cd[g,]+1)), 
            log10(AKI_cortex.cd[g,]+1)/max(log10(AKI_cd[g,]+1)), 
            log10(CTRL_cortex.cd[g,]+1)/max(log10(ctrl_cd[g,]+1)), 
            log10(IRL_cortex.cd[g,]+1)/max(log10(irl_cd[g,]+1)),
            log10(CIS_interface.cd[g,]+1)/max(log10(CIS_cd[g,]+1)), 
            log10(AKI_interface.cd[g,]+1)/max(log10(AKI_cd[g,]+1)), 
            log10(CTRL_interface.cd[g,]+1)/max(log10(ctrl_cd[g,]+1)), 
            log10(IRL_interface.cd[g,]+1)/max(log10(irl_cd[g,]+1)),
            log10( CIS_medulla.cd[g,]+1)/max(log10(CIS_cd[g,]+1)), 
            log10(AKI_medulla.cd[g,]+1)/max(log10(AKI_cd[g,]+1)), 
            log10(CTRL_medulla.cd[g,]+1)/max(log10(ctrl_cd[g,]+1)), 
            log10(IRL_medulla.cd[g,]+1)/max(log10(irl_cd[g,]+1)) ),
  
 group = factor(c(rep("CIS_Cortex", length(CIS_cortex.cd[g,])),
                 rep("AKI_Cortex", length(AKI_cortex.cd[g,])),
                 rep("CTRL_Cortex", length(CTRL_cortex.cd[g,])),
                 rep("IRL_Cortex", length(IRL_cortex.cd[g,])),
                 rep("CIS_Interface", length(CIS_interface.cd[g,])),
                 rep("AKI_Interface", length(AKI_interface.cd[g,])),
                 rep("CTRL_Interface", length(CTRL_interface.cd[g,])),
                 rep("IRL_Interface", length(IRL_interface.cd[g,])),
                 rep("CIS_Medulla", length(CIS_medulla.cd[g,])),
                 rep("AKI_Medulla", length(AKI_medulla.cd[g,])),
                 rep("CTRL_Medulla", length(CTRL_medulla.cd[g,])),
                 rep("IRL_Medulla", length(IRL_medulla.cd[g,])) ))
 )

# data <- data.frame(
#   value = c(CIS_cortex.cd[g,], AKI_cortex.cd[g,], CTRL_cortex.cd[g,], IRL_cortex.cd[g,],
#             CIS_interface.cd[g,], AKI_interface.cd[g,], CTRL_interface.cd[g,], IRL_interface.cd[g,],
#             CIS_medulla.cd[g,], AKI_medulla.cd[g,], CTRL_medulla.cd[g,], IRL_medulla.cd[g,]),
#   group = factor(c(rep("CIS_Cortex", length(CIS_cortex.cd[g,])),
#                    rep("AKI_Cortex", length(AKI_cortex.cd[g,])),
#                    rep("CTRL_Cortex", length(CTRL_cortex.cd[g,])),
#                    rep("IRL_Cortex", length(IRL_cortex.cd[g,])),
#                    rep("CIS_Interface", length(CIS_interface.cd[g,])),
#                    rep("AKI_Interface", length(AKI_interface.cd[g,])),
#                    rep("CTRL_Interface", length(CTRL_interface.cd[g,])),
#                    rep("IRL_Interface", length(IRL_interface.cd[g,])),
#                    rep("CIS_Medulla", length(CIS_medulla.cd[g,])),
#                    rep("AKI_Medulla", length(AKI_medulla.cd[g,])),
#                    rep("CTRL_Medulla", length(CTRL_medulla.cd[g,])),
#                    rep("IRL_Medulla", length(IRL_medulla.cd[g,]))
#                    
#   ) ) )


data$group<-factor(data$group, levels= c("CIS_Cortex", "AKI_Cortex", "CTRL_Cortex", "IRL_Cortex",
                                         "CIS_Interface", "AKI_Interface", "CTRL_Interface", "IRL_Interface",
                                         "CIS_Medulla", "AKI_Medulla", "CTRL_Medulla", "IRL_Medulla"))

# # Load ggplot2
library(ggplot2)

# Plot the violin plot
ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7, width=1) +  # Violin plot with transparency
  geom_boxplot(width = 0.05, outlier.shape = NA,   # Boxplot with reduced width
               alpha = 0.9, color = "black") +    # Boxplot with black outline

  theme_minimal() +          # Use a clean theme
  labs(x = "Group", y = "Norm. log(CPM)", title= g) +
  scale_fill_brewer(palette = "Set3") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +# Optional color palette
  theme(legend.position = "none")

# Also make a bubbule plot for the gene of interest
ggplot(data, aes(x = group, y = 'Slc34a1', size = value)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(2, 10)) + # Adjust bubble size range
  scale_color_gradient(low = "blue", high = "red") + # Adjust color gradient
  labs(
    title = "Bubble Plot of Gene Expression",
    x = "Sample",
    y = "Gene",
    size = "Expression Level",
    color = "Expression Level"
  ) +
  theme_minimal()

########################## PLOT 3 ##############################################
#Gene expression matrix from the all the kidney datasets
cd<-cbind(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),],
          CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),],
          CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),],
          CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),],
          AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl1$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl2$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl3$gexp[unique(rownames(CIS_0h$gexp)),],
          ctrl4$gexp[unique(rownames(CIS_0h$gexp)),],
          irl1$gexp[unique(rownames(CIS_0h$gexp)),],
          irl2$gexp[unique(rownames(CIS_0h$gexp)),],
          irl3$gexp[unique(rownames(CIS_0h$gexp)),],
          irl4$gexp[unique(rownames(CIS_0h$gexp)),] )

#CPM Normalization
cd<-MERINGUE::normalizeCounts(cd,log=FALSE)

#Subsetting the count matrix into different compartments:
cd.cortex<-cd[,all_spots$Tissue=='Cortex']
cd.interface<-cd[,all_spots$Tissue=='Interface']
cd.medulla<-cd[,all_spots$Tissue=='Medulla']
cd.other<-cd[,all_spots$Tissue=='Other']

# Plot the violin plot
g<-'Slc34a1'
data <- data.frame(
  value = c(cd.cortex[g,], cd.interface[g,], cd.medulla[g,], cd.other[g,]),
  group = factor(c(rep("Cortex", length(cd.cortex[g,])),
                   rep("Interface", length(cd.interface[g,])),
                   rep("Medulla", length(cd.medulla[g,])),
                   rep("Other", length(cd.other[g,]))
  ) ) )

library(ggplot2)
ggplot(data, aes(x = group, y = log(value)/max(log(value)), fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7, width=1) +  # Violin plot with transparency
  geom_boxplot(width = 0.05, outlier.shape = NA,   # Boxplot with reduced width
               alpha = 0.9, color = "black") +    # Boxplot with black outline
  
  theme_minimal() +          # Use a clean theme
  labs(x = "Group", y = "Value") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_fill_brewer(palette = "Set3") + # Optional color palette
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  theme(legend.position = "none")

##xxxx

g<-'Cyp7b1'
data <- data.frame(
  value = c(cd.cortex[g,], cd.interface[g,], cd.medulla[g,], cd.other[g,]),
  group = factor(c(rep("Cortex", length(cd.cortex[g,])),
                   rep("Interface", length(cd.interface[g,])),
                   rep("Medulla", length(cd.medulla[g,])),
                   rep("Other", length(cd.other[g,]))
  ) ) )

library(ggplot2)
ggplot(data, aes(x = group, y = log(value)/max(log(value)), fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7, width=1) +  # Violin plot with transparency
  geom_boxplot(width = 0.05, outlier.shape = NA,   # Boxplot with reduced width
               alpha = 0.9, color = "black") +    # Boxplot with black outline
  
  theme_minimal() +          # Use a clean theme
  labs(x = "Group", y = "Value") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_fill_brewer(palette = "Set3") + # Optional color palette
  theme(legend.position = "none")

g<-'Slc14a2'
data <- data.frame(
  value = c(cd.cortex[g,], cd.interface[g,], cd.medulla[g,]),
  group = factor(c(rep("Cortex", length(cd.cortex[g,])),
                   rep("Interface", length(cd.interface[g,])),
                   rep("Medulla", length(cd.medulla[g,]))
  ) ) )


range(log10(data$value+1))
hist(data$value)

library(ggplot2)
ggplot(data, aes(x = group, y = log(value+1), fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7, width=1) +  # Violin plot with transparency
  geom_boxplot(width = 0.05, outlier.shape = NA,   # Boxplot with reduced width
               alpha = 0.9, color = "black") +    # Boxplot with black outline
  
  theme_minimal() +          # Use a clean theme
  labs(y ="log(CPM)", title = g) + 
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_fill_brewer(palette = "Set3") + # Optional color palette
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  theme(legend.position = "none")
