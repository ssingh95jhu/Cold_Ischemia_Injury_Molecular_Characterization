#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

library(openxlsx)
library(ggplot2)
#library(dplyr)
#library(reshape2)

####### COMPARTMENTAL ENRICHED PATHWAYS BARPLOTS ###############################
#Note: Compartmental enriched pathways are ordered by their normalized enrichment
#scores (NES) obtained from gene set enrichemnt analysis.
#CIS Cortex
CIS_Cortex<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Cortex")
CIS_Cortex$Description<-gsub(" - Mus.*","",CIS_Cortex$Description)
df_Cortex<-data.frame(Pathways=CIS_Cortex$Description, NES=CIS_Cortex$NES, Enrichment="Positive")
df_Cortex[which(df_Cortex$NES<0),]$Enrichment<-"Negative"
#rownames(df_Cortex)<-CIS_Cortex$Description  
head(df_Cortex)

df_Cortex$Pathways=factor(df_Cortex$Pathways, levels = rev(df_Cortex$Pathways))

pdf("Figures/FIgure3/pdfs/GSEA_NES_Cortex_Barplots.pdf", width=8, height=16)
ggplot(df_Cortex, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (KEGG)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Cortex)") +
  theme_minimal() + scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
    dev.off()

#CIS Interface (Interface is same as Outer Medulla)
CIS_Interface<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Interface")
CIS_Interface$Description<-gsub(" - Mus.*","",CIS_Interface$Description)
df_Interface<-data.frame(Pathways=CIS_Interface$Description, NES=CIS_Interface$NES, Enrichment="Positive")
df_Interface[which(df_Interface$NES<0),]$Enrichment<-"Negative"
#rownames(df_Interface)<-CIS_Interface$Description  
head(df_Interface)

df_Interface$Pathways=factor(df_Interface$Pathways, levels = rev(df_Interface$Pathways))

pdf("Figures/FIgure3/pdfs/GSEA_NES_Interface_Barplots.pdf", width=8, height=16)
ggplot(df_Interface, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (KEGG)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Interface)") +
  theme_minimal() + scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"), # Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()

#CIS Medulla (Medulla is the same as Inner Medulla)
CIS_Medulla<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Medulla")
CIS_Medulla$Description<-gsub(" - Mus.*","",CIS_Medulla$Description)
df_Medulla<-data.frame(Pathways=CIS_Medulla$Description, NES=CIS_Medulla$NES, Enrichment="Positive")
df_Medulla[which(df_Medulla$NES<0),]$Enrichment<-"Negative"
#rownames(df_Medulla)<-CIS_Medulla$Description  
head(df_Medulla)

df_Medulla$Pathways=factor(df_Medulla$Pathways, levels = rev(df_Medulla$Pathways))

pdf("Figures/FIgure3/pdfs/GSEA_NES_Medulla_Barplots.pdf", width=8, height=16)
ggplot(df_Medulla, aes(x = Pathways, y = NES, fill=Enrichment)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Positive" = "orangered", "Negative" = "deepskyblue")) + 
  labs(x = "Enrcihed Pathways (KEGG)", y = "Normalized Enrichment Score", title = "Gene Set Enrichment Analysis (Medulla)") +
  theme_minimal() + scale_y_continuous(limits=c(-7,7), breaks = c(seq(-7, 7, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"),# Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()
