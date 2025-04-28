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

####### ALL AKI PATHWAYS (KEGG) ###########################################################
#Cortex
AKI_Cortex<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_KEGG_Pathways.xlsx", sheet="AKI_Cortex")
AKI_Cortex$Description<-gsub(" - Mus.*","",AKI_Cortex$Description)
df_Cortex<-data.frame(Pathways=AKI_Cortex$Description, NES=AKI_Cortex$NES, Enrichment="Positive")
df_Cortex[which(df_Cortex$NES<0),]$Enrichment<-"Negative"
#rownames(df_Cortex)<-AKI_Cortex$Description  
head(df_Cortex)


df_Cortex$Pathways=factor(df_Cortex$Pathways, levels = rev(df_Cortex$Pathways))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_NES_Cortex_Barplots.pdf", width=8, height=16)
ggplot(df_Cortex, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (KEGG)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Cortex)") +
  theme_minimal() + scale_y_continuous(limits=c(-10,10), breaks = c(seq(-10, 10, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 9, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
    dev.off()

#Interface
AKI_Interface<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_KEGG_Pathways.xlsx", sheet="AKI_Interface")
AKI_Interface$Description<-gsub(" - Mus.*","",AKI_Interface$Description)
df_Interface<-data.frame(Pathways=AKI_Interface$Description, NES=AKI_Interface$NES, Enrichment="Positive")
df_Interface[which(df_Interface$NES<0),]$Enrichment<-"Negative"
#rownames(df_Interface)<-AKI_Interface$Description  
head(df_Interface)

df_Interface$Pathways=factor(df_Interface$Pathways, levels = rev(df_Interface$Pathways))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_NES_Interface_Barplots.pdf", width=8, height=16)
ggplot(df_Interface, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (KEGG)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Interface)") +
  theme_minimal() + scale_y_continuous(limits=c(-10,10), breaks = c(seq(-10, 10, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 9, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()

#Medulla
AKI_Medulla<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_KEGG_Pathways.xlsx", sheet="AKI_Medulla")
AKI_Medulla$Description<-gsub(" - Mus.*","",AKI_Medulla$Description)
df_Medulla<-data.frame(Pathways=AKI_Medulla$Description, NES=AKI_Medulla$NES, Enrichment="Positive")
df_Medulla[which(df_Medulla$NES<0),]$Enrichment<-"Negative"
#rownames(df_Medulla)<-AKI_Medulla$Description  
head(df_Medulla)

df_Medulla$Pathways=factor(df_Medulla$Pathways, levels = rev(df_Medulla$Pathways))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_NES_Medulla_Barplots.pdf", width=8, height=16)
ggplot(df_Medulla, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (KEGG)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Medulla)") +
  theme_minimal() + scale_y_continuous(limits=c(-10,10), breaks = c(seq(-10, 10, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 9, color="black"),# Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()

####### ALL AKI PATHWAYS (HALLAMRK) ###########################################################
#Cortex
AKI_Cortex<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_HALLMARK_Pathways.xlsx", sheet="AKI_Cortex")
AKI_Cortex$Description<-gsub("HALLMARK_","",AKI_Cortex$Description)
df_Cortex<-data.frame(Pathways=AKI_Cortex$Description, NES=AKI_Cortex$NES, Enrichment="Positive")
df_Cortex[which(df_Cortex$NES<0),]$Enrichment<-"Negative"
#rownames(df_Cortex)<-AKI_Cortex$Description  
head(df_Cortex)

df_Cortex$Pathways=factor(df_Cortex$Pathways, levels = rev(df_Cortex$Pathways))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_NES_Cortex_HALLAMARK_Barplots.pdf", width=8, height=16)
ggplot(df_Cortex, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (HALLMARK)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (AKI Cortex)") +
  theme_minimal() + scale_y_continuous(limits=c(-10,10), breaks = c(seq(-10, 10, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()

#Interface
AKI_Interface<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_HALLMARK_Pathways.xlsx", sheet="AKI_Interface")
AKI_Interface$Description<-gsub("HALLMARK_","",AKI_Interface$Description)
df_Interface<-data.frame(Pathways=AKI_Interface$Description, NES=AKI_Interface$NES, Enrichment="Positive")
df_Interface[which(df_Interface$NES<0),]$Enrichment<-"Negative"
#rownames(df_Interface)<-AKI_Interface$Description  
head(df_Interface)

df_Interface$Pathways=factor(df_Interface$Pathways, levels = rev(df_Interface$Pathways))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_NES_Interface_HALLMARK_Barplots.pdf", width=8, height=16)
ggplot(df_Interface, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (HALLMARK)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Interface)") +
  theme_minimal() + scale_y_continuous(limits=c(-10,10), breaks = c(seq(-10, 10, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()

#Medulla
AKI_Medulla<-read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_GSEA_HALLMARK_Pathways.xlsx", sheet="AKI_Medulla")
AKI_Medulla$Description<-gsub("HALLMARK_","",AKI_Medulla$Description)
df_Medulla<-data.frame(Pathways=AKI_Medulla$Description, NES=AKI_Medulla$NES, Enrichment="Positive")
df_Medulla[which(df_Medulla$NES<0),]$Enrichment<-"Negative"
#rownames(df_Medulla)<-AKI_Medulla$Description  
head(df_Medulla)

df_Medulla$Pathways=factor(df_Medulla$Pathways, levels = rev(df_Medulla$Pathways))

pdf("Figures/Figure4S/pdfs/GSEA_AKI_NES_Medulla_HALLAMRK_Barplots.pdf", width=8, height=16)
ggplot(df_Medulla, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (KEGG)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Medulla)") +
  theme_minimal() + scale_y_continuous(limits=c(-10,10), breaks = c(seq(-10, 10, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"),# Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()

