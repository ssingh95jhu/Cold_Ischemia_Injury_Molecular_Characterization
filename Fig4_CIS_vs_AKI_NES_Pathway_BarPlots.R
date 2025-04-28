#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(cowplot)
library(viridis) #color blind friendly
library(RColorBrewer) #for RColorbrewer
library(reshape2)

#################### CREATE WORKBOOKS ##########################################
wb1<-createWorkbook()

######################## 1. PATHWAYS VISUALIZATION ################################
#######Extracting all saved pathways
### 1.CORTEX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cis_cortex_HALL<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_Hallmark_Pathways.xlsx", sheet = "CIS_Cortex")
rownames(cis_cortex_HALL)<-cis_cortex_HALL$ID
head(cis_cortex_HALL)

aki_cortex_HALL<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_Hallmark_Pathways.xlsx", sheet = "AKI_Cortex")
rownames(aki_cortex_HALL)<-aki_cortex_HALL$ID
head(aki_cortex_HALL)

cis_cortex_KEGG<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet = "CIS_Cortex")
rownames(cis_cortex_KEGG)<-cis_cortex_KEGG$ID
head(cis_cortex_KEGG)

aki_cortex_KEGG<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_KEGG_Pathways.xlsx", sheet = "AKI_Cortex")
rownames(aki_cortex_KEGG)<-aki_cortex_KEGG$ID
head(aki_cortex_KEGG)


### 2.INTERFACE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cis_interface_HALL<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_Hallmark_Pathways.xlsx", sheet = "CIS_Interface")
rownames(cis_interface_HALL)<-cis_interface_HALL$ID
head(cis_interface_HALL)

aki_interface_HALL<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_Hallmark_Pathways.xlsx", sheet = "AKI_Interface")
rownames(aki_interface_HALL)<-aki_interface_HALL$ID
head(aki_interface_HALL)

cis_interface_KEGG<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet = "CIS_Interface")
rownames(cis_interface_KEGG)<-cis_interface_KEGG$ID
head(cis_interface_KEGG)

aki_interface_KEGG<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_KEGG_Pathways.xlsx", sheet = "AKI_Interface")
rownames(aki_interface_KEGG)<-aki_interface_KEGG$ID
head(aki_interface_KEGG)


### 3. MEDULLA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cis_medulla_HALL<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_Hallmark_Pathways.xlsx", sheet = "CIS_Medulla")
rownames(cis_medulla_HALL)<-cis_medulla_HALL$ID
head(cis_medulla_HALL)

aki_medulla_HALL<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_Hallmark_Pathways.xlsx", sheet = "AKI_Medulla")
rownames(aki_medulla_HALL)<-aki_medulla_HALL$ID
head(aki_medulla_HALL)

cis_medulla_KEGG<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet = "CIS_Medulla")
rownames(cis_medulla_KEGG)<-cis_medulla_KEGG$ID
head(cis_medulla_KEGG)

aki_medulla_KEGG<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_KEGG_Pathways.xlsx", sheet = "AKI_Medulla")
rownames(aki_medulla_KEGG)<-aki_medulla_KEGG$ID
head(aki_medulla_KEGG)

## CIS CORTEX vs AKI CORTEX ( COAVRYING KEGG) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_KEGG<-intersect(rownames(cis_cortex_KEGG), rownames(aki_cortex_KEGG) )
df<-data.frame(cis_cortex_KEGG[comm_KEGG,]$NES, aki_cortex_KEGG[comm_KEGG,]$NES)
df$des<-cis_cortex_KEGG[comm_KEGG,]$Description
df$Sdes<-gsub(" - Mus musculus.*", "", df$des)

colnames(df)<-c('CIS_cortex_KEGG_NES','AKI_cortex_KEGG_NES', 'CIS_AKI_KEGG_Path', 'ShortPath')
head(df)


df_covary<-df[which((df$CIS_cortex_KEGG_NES*df$AKI_cortex_KEGG_NES)>0),]
head(df_covary)

CIS_Cor_KEGG<-data.frame(KEGG_Path=df_covary$ShortPath, CIS_Cortex_NES=df_covary$CIS_cortex_KEGG_NES, Enrichment="positive")
#CIS_Cor_KEGG_plot<-rbind(head(CIS_Cor_KEGG), tail(CIS_Cor_KEGG))
CIS_Cor_KEGG_up<-CIS_Cor_KEGG[which(CIS_Cor_KEGG$CIS_Cortex_NES>0),]
CIS_Cor_KEGG_down<-CIS_Cor_KEGG[which(CIS_Cor_KEGG$CIS_Cortex_NES<0),]
CIS_Cor_KEGG_plot<-rbind(CIS_Cor_KEGG_up, CIS_Cor_KEGG_down)
CIS_Cor_KEGG_plot$Enrichment<-ifelse(CIS_Cor_KEGG_plot$CIS_Cortex_NES<0, "negative", "positive")
CIS_Cor_KEGG_plot$CIS_Cortex_NES<-ifelse(CIS_Cor_KEGG_plot$CIS_Cortex_NES<0, -CIS_Cor_KEGG_plot$CIS_Cortex_NES, CIS_Cor_KEGG_plot$CIS_Cortex_NES)

AKI_Cor_KEGG<-data.frame(KEGG_Path=df_covary$ShortPath, AKI_Cortex_NES=df_covary$AKI_cortex_KEGG_NES, Enrichment="positive")
# AKI_Cor_KEGG_plot<-rbind(head(AKI_Cor_KEGG), tail(AKI_Cor_KEGG))
AKI_Cor_KEGG_up<-AKI_Cor_KEGG[which(AKI_Cor_KEGG$AKI_Cortex_NES>0),]
AKI_Cor_KEGG_down<-AKI_Cor_KEGG[which(AKI_Cor_KEGG$AKI_Cortex_NES<0),]
AKI_Cor_KEGG_plot<-rbind(AKI_Cor_KEGG_up, AKI_Cor_KEGG_down)
AKI_Cor_KEGG_plot$Enrichment<-ifelse(AKI_Cor_KEGG_plot$AKI_Cortex_NES<0, "negative", "positive")
AKI_Cor_KEGG_plot$AKI_Cortex_NES<-ifelse(AKI_Cor_KEGG_plot$AKI_Cortex_NES>0, -AKI_Cor_KEGG_plot$AKI_Cortex_NES, AKI_Cor_KEGG_plot$AKI_Cortex_NES)

colnames(CIS_Cor_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
colnames(AKI_Cor_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
CIS_AKI_Cortex_KEGG_plot<-rbind(CIS_Cor_KEGG_plot,AKI_Cor_KEGG_plot)

pdf("Figures/Figure4/pdfs/CIS_vs_AKI_GSEA_NES_Cortex_Barplots.pdf", width=8, height=16)
ggplot(CIS_AKI_Cortex_KEGG_plot, aes(x = rev(factor(KEGG_Path, levels=unique(KEGG_Path))), y = NES , fill = Enrichment)) +
  #geom_bar(stat = "identity", position = "identity", width=0.7, alpha =0.5) + coord_flip() + 
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = rev(unique(CIS_AKI_Cortex_KEGG_plot$KEGG_Path))) +
  #geom_text(aes(label = KEGG_Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
  labs(x="",y = "Normalized Enrichment Score", title = "CIS Cortex vs AKI Cortex") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) 
dev.off()

#Disjoint Pathways ************************************************************
df_disjoint<-df[which((df$CIS_cortex_KEGG_NES*df$AKI_cortex_KEGG_NES)<0),]
head(df_disjoint)

CIS_Cor_KEGG<-data.frame(KEGG_Path=df_disjoint$ShortPath, CIS_Cortex_NES=df_disjoint$CIS_cortex_KEGG_NES, Enrichment="positive")
#CIS_Cor_KEGG_plot<-rbind(head(CIS_Cor_KEGG), tail(CIS_Cor_KEGG))
CIS_Cor_KEGG_up<-CIS_Cor_KEGG[which(CIS_Cor_KEGG$CIS_Cortex_NES>0),]
CIS_Cor_KEGG_down<-CIS_Cor_KEGG[which(CIS_Cor_KEGG$CIS_Cortex_NES<0),]
CIS_Cor_KEGG_plot<-rbind(CIS_Cor_KEGG_up, CIS_Cor_KEGG_down)
CIS_Cor_KEGG_plot$Enrichment<-ifelse(CIS_Cor_KEGG_plot$CIS_Cortex_NES<0, "negative", "positive")
CIS_Cor_KEGG_plot$CIS_Cortex_NES<-ifelse(CIS_Cor_KEGG_plot$CIS_Cortex_NES<0, -CIS_Cor_KEGG_plot$CIS_Cortex_NES, CIS_Cor_KEGG_plot$CIS_Cortex_NES)

AKI_Cor_KEGG<-data.frame(KEGG_Path=df_disjoint$ShortPath, AKI_Cortex_NES=df_disjoint$AKI_cortex_KEGG_NES, Enrichment="positive")
#AKI_Cor_KEGG_plot<-rbind(head(AKI_Cor_KEGG), tail(AKI_Cor_KEGG))
AKI_Cor_KEGG_up<-AKI_Cor_KEGG[which(AKI_Cor_KEGG$AKI_Cortex_NES>0),]
AKI_Cor_KEGG_down<-AKI_Cor_KEGG[which(AKI_Cor_KEGG$AKI_Cortex_NES<0),]
AKI_Cor_KEGG_plot<-rbind(AKI_Cor_KEGG_up, AKI_Cor_KEGG_down)
AKI_Cor_KEGG_plot$Enrichment<-ifelse(AKI_Cor_KEGG_plot$AKI_Cortex_NES<0, "negative", "positive")
AKI_Cor_KEGG_plot$AKI_Cortex_NES<-ifelse(AKI_Cor_KEGG_plot$AKI_Cortex_NES>0, -AKI_Cor_KEGG_plot$AKI_Cortex_NES, AKI_Cor_KEGG_plot$AKI_Cortex_NES)

colnames(CIS_Cor_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
colnames(AKI_Cor_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
CIS_AKI_Cortex_KEGG_plot<-rbind(CIS_Cor_KEGG_plot,AKI_Cor_KEGG_plot)
#CIS_AKI_Cortex_KEGG_plot<-rbind(AKI_Cor_KEGG_plot,CIS_Cor_KEGG_plot)

pdf("Figures/Figure4S/pdfs/CIS_vs_AKI_GSEA_NES_Cortex_Barplots.pdf", width=8, height=16)
ggplot(CIS_AKI_Cortex_KEGG_plot, aes(x = factor(KEGG_Path, levels=unique(KEGG_Path)), y = NES , fill = Enrichment)) +
  #geom_bar(stat = "identity", position = "identity", width=0.7, alpha =0.5) + coord_flip() + 
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = unique(CIS_AKI_Cortex_KEGG_plot$KEGG_Path)) +
  #geom_text(aes(label = KEGG_Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
  labs(y = "Normalized Enrichment Score", x="KEGG_Pathways", title = "CISvsAKI_Cortex_Disjoint") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()


# CIS INTERFACE vs AKI INTERFACE (KEGG) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_KEGG<-intersect(rownames(cis_interface_KEGG), rownames(aki_interface_KEGG) )
df<-data.frame(cis_interface_KEGG[comm_KEGG,]$NES, aki_interface_KEGG[comm_KEGG,]$NES)
df$des<-cis_interface_KEGG[comm_KEGG,]$Description
df$Sdes<-gsub(" - Mus musculus.*", "", df$des)

colnames(df)<-c('CIS_interface_KEGG_NES','AKI_interface_KEGG_NES', 'CIS_AKI_KEGG_Path', 'ShortPath')
head(df)

##Covarying Pathways
df_covary<-df[which((df$CIS_interface_KEGG_NES*df$AKI_interface_KEGG_NES)>0),]
head(df_covary)

CIS_Int_KEGG<-data.frame(KEGG_Path=df_covary$ShortPath, CIS_Interface_NES=df_covary$CIS_interface_KEGG_NES, Enrichment="positive")
#CIS_Med_KEGG_plot<-rbind(head(CIS_Med_KEGG), tail(CIS_Med_KEGG))
CIS_Int_KEGG_up<-CIS_Int_KEGG[which(CIS_Int_KEGG$CIS_Interface_NES>0),]
CIS_Int_KEGG_down<-CIS_Int_KEGG[which(CIS_Int_KEGG$CIS_Interface_NES<0),]
CIS_Int_KEGG_plot<-rbind(CIS_Int_KEGG_up, CIS_Int_KEGG_down)
CIS_Int_KEGG_plot$Enrichment<-ifelse(CIS_Int_KEGG_plot$CIS_Interface_NES<0, "negative", "positive")
CIS_Int_KEGG_plot$CIS_Interface_NES<-ifelse(CIS_Int_KEGG_plot$CIS_Interface_NES<0, -CIS_Int_KEGG_plot$CIS_Interface_NES, CIS_Int_KEGG_plot$CIS_Interface_NES)

AKI_Int_KEGG<-data.frame(KEGG_Path=df_covary$ShortPath, AKI_Interface_NES=df_covary$AKI_interface_KEGG_NES, Enrichment="positive")
# AKI_Int_KEGG_plot<-rbind(head(AKI_Int_KEGG), tail(AKI_Int_KEGG))
AKI_Int_KEGG_up<-AKI_Int_KEGG[which(AKI_Int_KEGG$AKI_Interface_NES>0),]
AKI_Int_KEGG_down<-AKI_Int_KEGG[which(AKI_Int_KEGG$AKI_Interface_NES<0),]
AKI_Int_KEGG_plot<-rbind(AKI_Int_KEGG_up, AKI_Int_KEGG_down)
AKI_Int_KEGG_plot$Enrichment<-ifelse(AKI_Int_KEGG_plot$AKI_Interface_NES<0, "negative", "positive")
AKI_Int_KEGG_plot$AKI_Interface_NES<-ifelse(AKI_Int_KEGG_plot$AKI_Interface_NES>0, -AKI_Int_KEGG_plot$AKI_Interface_NES, AKI_Int_KEGG_plot$AKI_Interface_NES)

colnames(CIS_Int_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
colnames(AKI_Int_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
CIS_AKI_Interface_KEGG_plot<-rbind(CIS_Int_KEGG_plot,AKI_Int_KEGG_plot)

pdf("Figures/Figure4/pdfs/CIS_vs_AKI_GSEA_NES_Interface_Barplots.pdf", width=7, height=16)
ggplot(CIS_AKI_Interface_KEGG_plot, aes(x = rev(factor(KEGG_Path, levels=unique(KEGG_Path))), y = NES , fill = Enrichment)) +
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = rev(unique(CIS_AKI_Interface_KEGG_plot$KEGG_Path))) +
  #geom_text(aes(label = KEGG_Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
  labs(y = "Normalized Enrichment Score", x="", title = "CIS Outer Medulla vs AKI Outer Medulla") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position="none")
dev.off()

##Disjoint Pathways
#Disjoint Pathways ************************************************************
df_disjoint<-df[which((df$CIS_interface_KEGG_NES*df$AKI_interface_KEGG_NES)<0),]
head(df_disjoint)

CIS_Int_KEGG<-data.frame(KEGG_Path=df_disjoint$ShortPath, CIS_Interface_NES=df_disjoint$CIS_interface_KEGG_NES, Enrichment="positive")
#CIS_Int_KEGG_plot<-rbind(head(CIS_Int_KEGG), tail(CIS_Int_KEGG))
CIS_Int_KEGG_up<-CIS_Int_KEGG[which(CIS_Int_KEGG$CIS_Interface_NES>0),]
CIS_Int_KEGG_down<-CIS_Int_KEGG[which(CIS_Int_KEGG$CIS_Interface_NES<0),]
CIS_Int_KEGG_plot<-rbind(CIS_Int_KEGG_up, CIS_Int_KEGG_down)
CIS_Int_KEGG_plot$Enrichment<-ifelse(CIS_Int_KEGG_plot$CIS_Interface_NES<0, "negative", "positive")
CIS_Int_KEGG_plot$CIS_Interface_NES<-ifelse(CIS_Int_KEGG_plot$CIS_Interface_NES<0, -CIS_Int_KEGG_plot$CIS_Interface_NES, CIS_Int_KEGG_plot$CIS_Interface_NES)

AKI_Int_KEGG<-data.frame(KEGG_Path=df_disjoint$ShortPath, AKI_Interface_NES=df_disjoint$AKI_interface_KEGG_NES, Enrichment="positive")
#AKI_Int_KEGG_plot<-rbind(head(AKI_Int_KEGG), tail(AKI_Int_KEGG))
AKI_Int_KEGG_up<-AKI_Int_KEGG[which(AKI_Int_KEGG$AKI_Interface_NES>0),]
AKI_Int_KEGG_down<-AKI_Int_KEGG[which(AKI_Int_KEGG$AKI_Interface_NES<0),]
AKI_Int_KEGG_plot<-rbind(AKI_Int_KEGG_up, AKI_Int_KEGG_down)
AKI_Int_KEGG_plot$Enrichment<-ifelse(AKI_Int_KEGG_plot$AKI_Interface_NES<0, "negative", "positive")
AKI_Int_KEGG_plot$AKI_Interface_NES<-ifelse(AKI_Int_KEGG_plot$AKI_Interface_NES>0, -AKI_Int_KEGG_plot$AKI_Interface_NES, AKI_Int_KEGG_plot$AKI_Interface_NES)

colnames(CIS_Int_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
colnames(AKI_Int_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
CIS_AKI_Interface_KEGG_plot<-rbind(CIS_Int_KEGG_plot,AKI_Int_KEGG_plot)
#CIS_AKI_Interface_KEGG_plot<-rbind(AKI_Int_KEGG_plot,CIS_Int_KEGG_plot)

pdf("Figures/Figure4S/pdfs/CIS_vs_AKI_GSEA_NES_Interface_Barplots.pdf", width=8, height=16)
ggplot(CIS_AKI_Interface_KEGG_plot, aes(x = factor(KEGG_Path, levels=unique(KEGG_Path)), y = NES , fill = Enrichment)) +
  #geom_bar(stat = "identity", position = "identity", width=0.7, alpha =0.5) + coord_flip() + 
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = unique(CIS_AKI_Interface_KEGG_plot$KEGG_Path)) +
  #geom_text(aes(label = KEGG_Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-10,10), breaks = c(seq(-10, 10, by = 2), 0)) +   
  labs(y = "Normalized Enrichment Score", x="KEGG_Pathways", title = "CISvsAKI_Interface_Disjoint") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()



# CIS MEDULLA vs AKI MEDULLA (KEGG) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_KEGG<-intersect(rownames(cis_medulla_KEGG), rownames(aki_medulla_KEGG) )
df<-data.frame(cis_medulla_KEGG[comm_KEGG,]$NES, aki_medulla_KEGG[comm_KEGG,]$NES)
df$des<-cis_medulla_KEGG[comm_KEGG,]$Description
df$Sdes<-gsub(" - Mus musculus.*", "", df$des)

colnames(df)<-c('CIS_medulla_KEGG_NES','AKI_medulla_KEGG_NES', 'CIS_AKI_KEGG_Path', 'ShortPath')
head(df)

df_covary<-df[which((df$CIS_medulla_KEGG_NES*df$AKI_medulla_KEGG_NES)>0),]
head(df_covary)

CIS_Med_KEGG<-data.frame(KEGG_Path=df_covary$ShortPath, CIS_Medulla_NES=df_covary$CIS_medulla_KEGG_NES, Enrichment="positive")
#CIS_Med_KEGG_plot<-rbind(head(CIS_Med_KEGG), tail(CIS_Med_KEGG))
CIS_Med_KEGG_up<-CIS_Med_KEGG[which(CIS_Med_KEGG$CIS_Medulla_NES>0),]
CIS_Med_KEGG_down<-CIS_Med_KEGG[which(CIS_Med_KEGG$CIS_Medulla_NES<0),]
CIS_Med_KEGG_plot<-rbind(CIS_Med_KEGG_up, CIS_Med_KEGG_down)
CIS_Med_KEGG_plot$Enrichment<-ifelse(CIS_Med_KEGG_plot$CIS_Medulla_NES<0, "negative", "positive")
CIS_Med_KEGG_plot$CIS_Medulla_NES<-ifelse(CIS_Med_KEGG_plot$CIS_Medulla_NES<0, -CIS_Med_KEGG_plot$CIS_Medulla_NES, CIS_Med_KEGG_plot$CIS_Medulla_NES)

AKI_Med_KEGG<-data.frame(KEGG_Path=df_covary$ShortPath, AKI_Medulla_NES=df_covary$AKI_medulla_KEGG_NES, Enrichment="positive")
# AKI_Med_KEGG_plot<-rbind(head(AKI_Med_KEGG), tail(AKI_Med_KEGG))
AKI_Med_KEGG_up<-AKI_Med_KEGG[which(AKI_Med_KEGG$AKI_Medulla_NES>0),]
AKI_Med_KEGG_down<-AKI_Med_KEGG[which(AKI_Med_KEGG$AKI_Medulla_NES<0),]
AKI_Med_KEGG_plot<-rbind(AKI_Med_KEGG_up, AKI_Med_KEGG_down)
AKI_Med_KEGG_plot$Enrichment<-ifelse(AKI_Med_KEGG_plot$AKI_Medulla_NES<0, "negative", "positive")
AKI_Med_KEGG_plot$AKI_Medulla_NES<-ifelse(AKI_Med_KEGG_plot$AKI_Medulla_NES>0, -AKI_Med_KEGG_plot$AKI_Medulla_NES, AKI_Med_KEGG_plot$AKI_Medulla_NES)

colnames(CIS_Med_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
colnames(AKI_Med_KEGG_plot)<-c('KEGG_Path', 'NES', 'Enrichment')
CIS_AKI_Medulla_KEGG_plot<-rbind(CIS_Med_KEGG_plot,AKI_Med_KEGG_plot)

pdf("Figures/FIgure4/pdfs/CIS_vs_AKI_GSEA_NES_Medulla_Barplots.pdf", width=7, height=16)
ggplot(CIS_AKI_Medulla_KEGG_plot, aes(x = rev(factor(KEGG_Path, levels=unique(KEGG_Path))), y = NES , fill = Enrichment)) +
  geom_bar(stat = "identity", position = "identity") + coord_flip() + 
  scale_x_discrete(labels = rev(unique(CIS_AKI_Medulla_KEGG_plot$KEGG_Path))) +
  #geom_text(aes(label = KEGG_Path, y=NES*1.8),  color = "black", size = 3, fontface="bold") +
  scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) +   
  labs(y = "Normalized Enrichment Score", x="", title = "CIS Inner Medulla vs AKI Inner Medulla") +
  scale_fill_manual(values = c("negative" = "deepskyblue", "positive" = "orangered")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  #geom_vline(xintercept = 0.5, color = "black", linetype = "solid", size = 1) +
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position="none")
dev.off()