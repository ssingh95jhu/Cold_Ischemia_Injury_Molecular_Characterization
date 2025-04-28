#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

## load data
load("data/CIS_data.RData")
load("data/AKI_data.RData")
load("data/Rabb_ctrl.RData")

########################## ADDING LIBRARY ######################################.
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(UpSetR)
library(ggforce)

############ UPSET PLOT ##########
######################### Upset plot ###########################################
#Modifying fromList function to retain rownames as gene names:
#Courtesy: docmanny https://github.com/hms-dbmi/UpSetR/issues/85

fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

####################### CIS COKPARTMENT SPECIFIC UNIQUE DEGs ###################
wb1<-createWorkbook() #workbook for CIS KEGG 

########## CIS DEGs *******************************************************
cis_global_up<-read.xlsx("Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Global_Up")
rownames(cis_global_up)<-cis_global_up$Gene

cis_cortex_up<-read.xlsx("Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Cortex_Up")
rownames(cis_cortex_up)<-cis_cortex_up$Gene

cis_interface_up<-read.xlsx("Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Interface_Up")
rownames(cis_interface_up)<-cis_interface_up$Gene

cis_medulla_up<-read.xlsx("Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Medulla_Up")
rownames(cis_medulla_up)<-cis_medulla_up$Gene

cis_global_down<-read.xlsx("Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Global_Down")
rownames(cis_global_down)<-cis_global_down$Gene

cis_cortex_down<-read.xlsx("Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Cortex_Down")
rownames(cis_cortex_down)<-cis_cortex_down$Gene

cis_interface_down<-read.xlsx("Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Interface_Down")
rownames(cis_interface_down)<-cis_interface_down$Gene

cis_medulla_down<-read.xlsx("Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Medulla_Down")
rownames(cis_medulla_down)<-cis_medulla_down$Gene


#Upset plot to check overlap between pathways
#CIS UP GENES
cis.genes.up<-list(A=cis_cortex_up$Gene,
                  B=cis_interface_up$Gene,
                  C=cis_medulla_up$Gene)
names(cis.genes.up)<-c("CIScorU","CISintU","CISmedU")
upset(fromList(cis.genes.up) ,order.by = "freq")

lst<-fromList(cis.genes.up)

CIS_Ocortex_up<-lst[rownames(lst[rowSums(lst)==1 & lst$CIScorU==1,]),] 
head(rownames(CIS_Ocortex_up))

CIS_Ointerface_up<-lst[rownames(lst[rowSums(lst)==1 & lst$CISintU==1,]),]
head(rownames(CIS_Ointerface_up))

CIS_Omedulla_up<-lst[rownames(lst[rowSums(lst)==1 & lst$CISmedU==1,]),] 
head(rownames(CIS_Omedulla_up))

CIS_all_up<-lst[rownames(lst[rowSums(lst)==3,]),] 
head(rownames(CIS_all_up))

addWorksheet(wb1,"CIS_OCor_Up")
writeData(wb1, "CIS_OCor_Up", data.frame(Gene=rownames(CIS_Ocortex_up),
                                         Slope=cis_cortex_up[rownames(CIS_Ocortex_up),]$Slope))

addWorksheet(wb1,"CIS_OInt_Up")
writeData(wb1, "CIS_OInt_Up", data.frame(Gene=rownames(CIS_Ointerface_up),
                                         Slope=cis_interface_up[rownames(CIS_Ointerface_up),]$Slope))

addWorksheet(wb1,"CIS_OMed_Up")
writeData(wb1, "CIS_OMed_Up", data.frame(Gene=rownames(CIS_Omedulla_up),
                                         Slope=cis_medulla_up[rownames(CIS_Omedulla_up),]$Slope))

addWorksheet(wb1,"CIS_All_Up")
writeData(wb1, "CIS_All_Up", data.frame(Gene=rownames(CIS_all_up),
                                         Cor_Slope=cis_cortex_up[rownames(CIS_all_up),]$Slope,
                                         Int_Slope=cis_interface_up[rownames(CIS_all_up),]$Slope,
                                         Med_Slope=cis_medulla_up[rownames(CIS_all_up),]$Slope)  )

#saveWorkbook(wb1, file = "Tables/CIS_AKI_Upset_DEGs.xlsx", overwrite = TRUE)

#CIS DOWN GENES
cis.genes.down<-list(A=cis_cortex_down$Gene,
                   B=cis_interface_down$Gene,
                   C=cis_medulla_down$Gene)
names(cis.genes.down)<-c("CIScorD","CISintD","CISmedD")
upset(fromList(cis.genes.down) ,order.by = "freq")

lst<-fromList(cis.genes.down)

CIS_Ocortex_down<-lst[rownames(lst[rowSums(lst)==1 & lst$CIScorD==1,]),] 
head(rownames(CIS_Ocortex_down))

CIS_Ointerface_down<-lst[rownames(lst[rowSums(lst)==1 & lst$CISintD==1,]),]
head(rownames(CIS_Ointerface_down))

CIS_Omedulla_down<-lst[rownames(lst[rowSums(lst)==1 & lst$CISmedD==1,]),] 
head(rownames(CIS_Omedulla_down))

CIS_all_down<-lst[rownames(lst[rowSums(lst)==3,]),] 
head(rownames(CIS_all_down))


addWorksheet(wb1,"CIS_OCor_Down")
writeData(wb1, "CIS_OCor_Down", data.frame(Gene=rownames(CIS_Ocortex_down),
                                         Slope=cis_cortex_down[rownames(CIS_Ocortex_down),]$Slope))

addWorksheet(wb1,"CIS_OInt_Down")
writeData(wb1, "CIS_OInt_Down", data.frame(Gene=rownames(CIS_Ointerface_down),
                                         Slope=cis_interface_down[rownames(CIS_Ointerface_down),]$Slope))

addWorksheet(wb1,"CIS_OMed_Down")
writeData(wb1, "CIS_OMed_Down", data.frame(Gene=rownames(CIS_Omedulla_down),
                                         Slope=cis_medulla_down[rownames(CIS_Omedulla_down),]$Slope))

addWorksheet(wb1,"CIS_All_Down")
writeData(wb1, "CIS_All_Down", data.frame(Gene=rownames(CIS_all_down),
                                        Cor_Slope=cis_cortex_down[rownames(CIS_all_down),]$Slope,
                                        Int_Slope=cis_interface_down[rownames(CIS_all_down),]$Slope,
                                        Med_Slope=cis_medulla_down[rownames(CIS_all_down),]$Slope)  )


####################### AKI COKPARTMENT SPECIFIC UNIQUE DEGs ##############
aki_cortex_up<-read.xlsx("Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Cortex_Up")
rownames(aki_cortex_up)<-aki_cortex_up$Gene

aki_interface_up<-read.xlsx("Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Interface_Up")
rownames(aki_interface_up)<-aki_interface_up$Gene

aki_medulla_up<-read.xlsx("Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Medulla_Up")
rownames(aki_medulla_up)<-aki_medulla_up$Gene

aki_cortex_down<-read.xlsx("Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Cortex_Down")
rownames(aki_cortex_down)<-aki_cortex_down$Gene

aki_interface_down<-read.xlsx("Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Interface_Down")
rownames(aki_interface_down)<-aki_interface_down$Gene

aki_medulla_down<-read.xlsx("Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Medulla_Down")
rownames(aki_medulla_down)<-aki_medulla_down$Gene


#Upset plot to check overlap between pathways
#AKI UP GENES
aki.genes.up<-list(A=aki_cortex_up$Gene,
                   B=aki_interface_up$Gene,
                   C=aki_medulla_up$Gene)
names(aki.genes.up)<-c("AKIcorU","AKIintU","AKImedU")
upset(fromList(aki.genes.up) ,order.by = "freq")

lst<-fromList(aki.genes.up)

AKI_Ocortex_up<-lst[rownames(lst[rowSums(lst)==1 & lst$AKIcorU==1,]),] 
head(rownames(AKI_Ocortex_up))

AKI_Ointerface_up<-lst[rownames(lst[rowSums(lst)==1 & lst$AKIintU==1,]),]
head(rownames(AKI_Ointerface_up))

AKI_Omedulla_up<-lst[rownames(lst[rowSums(lst)==1 & lst$AKImedU==1,]),] 
head(rownames(AKI_Omedulla_up))

AKI_all_up<-lst[rownames(lst[rowSums(lst)==3,]),] 
head(rownames(AKI_all_up))

addWorksheet(wb1,"AKI_OCor_Up")
writeData(wb1, "AKI_OCor_Up", data.frame(Gene=rownames(AKI_Ocortex_up),
                                         Slope=aki_cortex_up[rownames(AKI_Ocortex_up),]$Slope))

addWorksheet(wb1,"AKI_OInt_Up")
writeData(wb1, "AKI_OInt_Up", data.frame(Gene=rownames(AKI_Ointerface_up),
                                         Slope=aki_interface_up[rownames(AKI_Ointerface_up),]$Slope))

addWorksheet(wb1,"AKI_OMed_Up")
writeData(wb1, "AKI_OMed_Up", data.frame(Gene=rownames(AKI_Omedulla_up),
                                         Slope=aki_medulla_up[rownames(AKI_Omedulla_up),]$Slope))

addWorksheet(wb1,"AKI_All_Up")
writeData(wb1, "AKI_All_Up", data.frame(Gene=rownames(AKI_all_up),
                                        Cor_Slope=aki_cortex_up[rownames(AKI_all_up),]$Slope,
                                        Int_Slope=aki_interface_up[rownames(AKI_all_up),]$Slope,
                                        Med_Slope=aki_medulla_up[rownames(AKI_all_up),]$Slope)  )

#AKI DOWN GENES
aki.genes.down<-list(A=aki_cortex_down$Gene,
                   B=aki_interface_down$Gene,
                   C=aki_medulla_down$Gene)
names(aki.genes.down)<-c("AKIcorD","AKIintD","AKImedD")
upset(fromList(aki.genes.down) ,order.by = "freq")

lst<-fromList(aki.genes.down)

AKI_Ocortex_down<-lst[rownames(lst[rowSums(lst)==1 & lst$AKIcorD==1,]),] 
head(rownames(AKI_Ocortex_down))

AKI_Ointerface_down<-lst[rownames(lst[rowSums(lst)==1 & lst$AKIintD==1,]),]
head(rownames(AKI_Ointerface_down))

AKI_Omedulla_down<-lst[rownames(lst[rowSums(lst)==1 & lst$AKImedD==1,]),] 
head(rownames(AKI_Omedulla_down))

AKI_all_down<-lst[rownames(lst[rowSums(lst)==3,]),] 
head(rownames(AKI_all_down))


addWorksheet(wb1,"AKI_OCor_Down")
writeData(wb1, "AKI_OCor_Down", data.frame(Gene=rownames(AKI_Ocortex_down),
                                         Slope=aki_cortex_down[rownames(AKI_Ocortex_down),]$Slope))

addWorksheet(wb1,"AKI_OInt_Down")
writeData(wb1, "AKI_OInt_Down", data.frame(Gene=rownames(AKI_Ointerface_down),
                                         Slope=aki_interface_down[rownames(AKI_Ointerface_down),]$Slope))

addWorksheet(wb1,"AKI_OMed_Down")
writeData(wb1, "AKI_OMed_Down", data.frame(Gene=rownames(AKI_Omedulla_down),
                                         Slope=aki_medulla_down[rownames(AKI_Omedulla_down),]$Slope))

addWorksheet(wb1,"AKI_All_Down")
writeData(wb1, "AKI_All_Down", data.frame(Gene=rownames(AKI_all_down),
                                        Cor_Slope=aki_cortex_down[rownames(AKI_all_down),]$Slope,
                                        Int_Slope=aki_interface_down[rownames(AKI_all_down),]$Slope,
                                        Med_Slope=aki_medulla_down[rownames(AKI_all_down),]$Slope)  )

saveWorkbook(wb1, file = "Tables/CIS_AKI_Upset_DEGs.xlsx", overwrite = TRUE)

################################################################################
cortex_sdmg<-read.xlsx("Tables/Male_vs_Female_CTRL_vs_AKI_Sham_DEGs.xlsx", sheet = "Cortex")
rownames(cortex_sdmg)<-cortex_sdmg[,1]
cortex_sdmg<-cortex_sdmg[which(abs(as.numeric(cortex_sdmg$log2FoldChange))>1),]

interface_sdmg<-read.xlsx("Tables/Male_vs_Female_CTRL_vs_AKI_Sham_DEGs.xlsx", sheet = "Interface")
rownames(interface_sdmg)<-interface_sdmg[,1]
interface_sdmg<-interface_sdmg[which(abs(as.numeric(interface_sdmg$log2FoldChange))>1),]

medulla_sdmg<-read.xlsx("Tables/Male_vs_Female_CTRL_vs_AKI_Sham_DEGs.xlsx", sheet = "Medulla")
rownames(medulla_sdmg)<-medulla_sdmg[,1]
medulla_sdmg<-medulla_sdmg[which(abs(as.numeric(medulla_sdmg$log2FoldChange))>1),]

a<-intersect(rownames(CIS_Ocortex_up),rownames(AKI_Ocortex_down))
b<-intersect(rownames(CIS_Ointerface_up),rownames(AKI_Ointerface_down))
c<-intersect(rownames(CIS_Omedulla_up),rownames(AKI_Omedulla_down))

table(a%in%cortex_sdmg)
table(b%in%interface_sdmg)
table(c%in%medulla_sdmg)


# cis_cortex_up[intersect(rownames(CIS_Ocortex_up),rownames(AKI_Ocortex_up)),]
# cis_interface_up[intersect(rownames(CIS_Ointerface_up),rownames(AKI_Ointerface_up)),]
# cis_medulla_up[intersect(rownames(CIS_Omedulla_up),rownames(AKI_Omedulla_up)),]
# 
cis_cortex_up[intersect(rownames(cis_cortex_up),rownames(aki_cortex_up)),]
# cis_interface_up[intersect(rownames(cis_interface_up),rownames(aki_interface_up)),]
# cis_medulla_up[intersect(rownames(cis_medulla_up),rownames(aki_medulla_up)),]
