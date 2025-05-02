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

####### ALL PATHWAYS ###########################################################
#Cortex
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

#Interface
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

#Medulla
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

###### ENERGY METABOLISM PATHWAYS ##############################################
#Cortex
CIS_Cortex<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Cortex")
CIS_Cortex$Description<-gsub(" - Mus.*","",CIS_Cortex$Description)

OXPHOS_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Oxidative phosphorylation"),]$NES
TCA_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Citrate cycle (TCA cycle)"),]$NES
Thermogenesis_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Thermogenesis"),]$NES

#Interface
CIS_Interface<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Interface")
CIS_Interface$Description<-gsub(" - Mus.*","",CIS_Interface$Description)

OXPHOS_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Oxidative phosphorylation"),]$NES
TCA_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Citrate cycle (TCA cycle)"),]$NES
Thermogenesis_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Thermogenesis"),]$NES

#Medulla
CIS_Medulla<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Medulla")
CIS_Medulla$Description<-gsub(" - Mus.*","",CIS_Medulla$Description)

OXPHOS_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Oxidative phosphorylation"),]$NES
TCA_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Citrate cycle (TCA cycle)"),]$NES
Thermogenesis_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Thermogenesis"),]$NES

data<-data.frame(Cortex=c(OXPHOS_Cor_NES, TCA_Cor_NES, Thermogenesis_Cor_NES),
                 Interface=c(OXPHOS_Int_NES, TCA_Int_NES, Thermogenesis_Int_NES),
                 Medulla=c(OXPHOS_Med_NES, TCA_Med_NES, Thermogenesis_Med_NES))
rownames(data)<-c("OXPHOS", "TCA", "Thermogenesis")

df_long <- melt(as.matrix(data), varnames = c("Pathway", "Region"), value.name = "NES")

ggplot(df_long, aes(x = Pathway, y = NES, fill = Region)) +
  geom_col(position = position_dodge2(width= 0.5, padding=0.1)) +  # Align bars properly
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_fill_manual(values = c("red", "blue", "green")) +# Avoid label overlap
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +  # Add horizontal line at y=0
  theme_minimal() +
  labs(title = "Energy Metabolism",
       y = "Normalized Enrichment Score",
       x = "Pathway") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color="black", size=12),
        axis.text.y = element_text(angle = 0, hjust = 0.5, color="black", size=12),
        plot.title = element_text(hjust = 0.5, size = 12, color="black" ), 
        axis.title.x = element_text(angle=0, hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5, angle = 90, vjust = 0.5, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks=element_line(color="black", linewidth=0.5),
        axis.ticks.length=unit(0.1,"cm"))  # Keep labels centered

###### INTERMEDIATE/SUBSTRATE GENERATING PATHWAYS ##############################
#Cortex
CIS_Cortex<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Cortex")
CIS_Cortex$Description<-gsub(" - Mus.*","",CIS_Cortex$Description)

FAO_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Fatty acid degradation"),]$NES
Tryptophan_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Tryptophan metabolism"),]$NES
Valine_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Valine, leucine and isoleucine degradation"),]$NES
Glycolysis_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Glycolysis / Gluconeogenesis"),]$NES
Glycolysis_Cor_NES<-0

#Interface
CIS_Interface<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Interface")
CIS_Interface$Description<-gsub(" - Mus.*","",CIS_Interface$Description)

FAO_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Fatty acid degradation"),]$NES
Tryptophan_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Tryptophan metabolism"),]$NES
Valine_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Valine, leucine and isoleucine degradation"),]$NES
Glycolysis_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Glycolysis / Gluconeogenesis"),]$NES
Glycolysis_Int_NES<-0

#Medulla
CIS_Medulla<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Medulla")
CIS_Medulla$Description<-gsub(" - Mus.*","",CIS_Medulla$Description)

FAO_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Fatty acid degradation"),]$NES
Tryptophan_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Tryptophan metabolism"),]$NES
Valine_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Valine, leucine and isoleucine degradation"),]$NES
Glycolysis_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Glycolysis / Gluconeogenesis"),]$NES

data<-data.frame(Cortex=c(FAO_Cor_NES, Tryptophan_Cor_NES, Valine_Cor_NES,Glycolysis_Cor_NES),
                 Interface=c(FAO_Int_NES, Tryptophan_Int_NES, Valine_Int_NES, Glycolysis_Int_NES),
                 Medulla=c(FAO_Med_NES, Tryptophan_Med_NES, Valine_Med_NES, Glycolysis_Med_NES))
rownames(data)<-c("FAO", "Tryptophan", "Valine", "Glycolysis")

df_long <- melt(as.matrix(data), varnames = c("Pathway", "Region"), value.name = "NES")

df_long[is.na(df_long$NES)] <- 0

ggplot(df_long, aes(x = Pathway, y = NES, fill = Region)) +
  geom_col(position = position_dodge2(width= 0.5, padding=0.1)) +  # Align bars properly
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_fill_manual(values = c("red", "blue", "green")) +# Avoid label overlap
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +  # Add horizontal line at y=0
  theme_minimal() +
  labs(title = "Substrate Metabolism",
       y = "Normalized Enrichment Score",
       x = "Pathway") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color="black", size=12),
        axis.text.y = element_text(angle = 0, hjust = 0.5, color="black", size=12),
        plot.title = element_text(hjust = 0.5, size = 12, color="black" ), 
        axis.title.x = element_text(angle=0, hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5, angle = 90, vjust = 0.5, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks=element_line(color="black", linewidth=0.5),
        axis.ticks.length=unit(0.1,"cm"))  # Keep labels centered

###### OXIDATIVE STRESS PATHWAYS ###############################################
#Cortex
CIS_Cortex<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Cortex")
CIS_Cortex$Description<-gsub(" - Mus.*","",CIS_Cortex$Description)

ROS_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Chemical carcinogenesis - reactive oxygen species"),]$NES
Nerve_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Pathways of neurodegeneration - multiple diseases"),]$NES

#Interface
CIS_Interface<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Interface")
CIS_Interface$Description<-gsub(" - Mus.*","",CIS_Interface$Description)

ROS_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Chemical carcinogenesis - reactive oxygen species"),]$NES
Nerve_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Pathways of neurodegeneration - multiple diseases"),]$NES

#Medulla
CIS_Medulla<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Medulla")
CIS_Medulla$Description<-gsub(" - Mus.*","",CIS_Medulla$Description)

ROS_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Chemical carcinogenesis - reactive oxygen species"),]$NES
Nerve_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Pathways of neurodegeneration - multiple diseases"),]$NES


data<-data.frame(Cortex=c(ROS_Cor_NES, Nerve_Cor_NES),
                 Interface=c(ROS_Int_NES, Nerve_Int_NES),
                 Medulla=c(ROS_Med_NES, Nerve_Med_NES))
rownames(data)<-c("ROS", "Nerve")

df_long <- melt(as.matrix(data), varnames = c("Pathway", "Region"), value.name = "NES")

ggplot(df_long, aes(x = Pathway, y = NES, fill = Region)) +
  geom_col(position = position_dodge2(width= 0.5, padding=0.1)) +  # Align bars properly
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_fill_manual(values = c("red", "blue", "green")) +# Avoid label overlap
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +  # Add horizontal line at y=0
  theme_minimal() +
  labs(title = "Oxidative Stress",
       y = "Normalized Enrichment Score",
       x = "Pathway") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color="black", size=12),
        axis.text.y = element_text(angle = 0, hjust = 0.5, color="black", size=12),
        plot.title = element_text(hjust = 0.5, size = 12, color="black" ), 
        axis.title.x = element_text(angle=0, hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5, angle = 90, vjust = 0.5, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks=element_line(color="black", linewidth=0.5),
        axis.ticks.length=unit(0.1,"cm"))  # Keep labels centered

###### PATHOLOGIC PATHWAYS ###############################################
#Cortex
CIS_Cortex<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Cortex")
CIS_Cortex$Description<-gsub(" - Mus.*","",CIS_Cortex$Description)

APC_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Antigen processing and presentation"),]$NES
Apoptosis_Cor_NES<-CIS_Cortex[which(CIS_Cortex$Description=="Apoptosis"),]$NES

#Interface
CIS_Interface<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Interface")
CIS_Interface$Description<-gsub(" - Mus.*","",CIS_Interface$Description)

APC_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Antigen processing and presentation"),]$NES
Apoptosis_Int_NES<-CIS_Interface[which(CIS_Interface$Description=="Apoptosis"),]$NES

#Medulla
CIS_Medulla<-read.xlsx("Supplementary_Tables/Limma_cpmCIS_GSEA_KEGG_Pathways.xlsx", sheet="CIS_Medulla")
CIS_Medulla$Description<-gsub(" - Mus.*","",CIS_Medulla$Description)

APC_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Antigen processing and presentation"),]$NES
Apoptosis_Med_NES<-CIS_Medulla[which(CIS_Medulla$Description=="Apoptosis"),]$NES

data<-data.frame(Cortex=c(APC_Cor_NES, Apoptosis_Cor_NES),
                 Interface=c(APC_Int_NES, Apoptosis_Int_NES),
                 Medulla=c(APC_Med_NES, Apoptosis_Med_NES))
rownames(data)<-c("APC", "Apoptosis")

df_long <- melt(as.matrix(data), varnames = c("Pathway", "Region"), value.name = "NES")

ggplot(df_long, aes(x = Pathway, y = NES, fill = Region)) +
  geom_col(position = position_dodge2(width= 0.5, padding=0.1)) +  # Align bars properly
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_fill_manual(values = c("red", "blue", "green")) +# Avoid label overlap
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +  # Add horizontal line at y=0
  theme_minimal() +
  labs(title = "Pathological Pathways",
       y = "Normalized Enrichment Score",
       x = "Pathway") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color="black", size=12),
        axis.text.y = element_text(angle = 0, hjust = 0.5, color="black", size=12),
        plot.title = element_text(hjust = 0.5, size = 12, color="black" ), 
        axis.title.x = element_text(angle=0, hjust = 0.5, size = 12),
        axis.title.y = element_text(hjust = 0.5, angle = 90, vjust = 0.5, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.ticks=element_line(color="black", linewidth=0.5),
        axis.ticks.length=unit(0.1,"cm"))  # Keep labels centered
