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

Imed_degs<-read.xlsx("Tables/Fig_4S_CTRL_vs_IRL_Compartmental_DEGs.xlsx", sheet=3)
rownames(Imed_degs)<-Imed_degs$Genes

goi<-c('Havcr1','Lcn2','Spp1','Hk1','Aldoa','Dlat','Pdha1','Pdhb','Got2','Ogdh',
       'Idh2','Idh3b','Mdh1','Mdh2','Sdha','Suclg2')
Imed_degs_sub<-Imed_degs[goi,]
df<-data.frame(Imed_degs_sub)

df$Trend<-'Up'
df[which(df$log2FoldChange<0),]$Trend<-'Down'
#df$Genes<-factor(df$Genes, levels=c('Havcr1','Lcn2','Spp1','Hk1','Aldoa','Dlat','Pdha1','Pdhb','Got2','Ogdh',
                                   # 'Idh2','Idh3b','Mdh1','Mdh2','Sdha','Suclg2'))

df$Genes<-factor(df$Genes, levels=rev(df$Genes))

pdf("Figures/FigureX/pdfs/CTRLvsAKI24_DEGs_InnerMedulla_Barplots.pdf", height=5, width=7)
ggplot(df, aes(x = Genes, y = as.numeric(log2FoldChange), fill=Trend)) + geom_bar(stat="identity") +
  scale_fill_manual(values = c("Up" = "orangered", "Down" = "deepskyblue")) + 
  labs(y = "Log2 Fold Change", title = "Native Kidney vs AKI24 Kidneys (Inner Medulla)") +
  theme_minimal() + scale_y_continuous(limits=c(-2,6), breaks = c(seq(-2, 6, by = 2), 0)) + 
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.5) + coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add box around plot
        axis.text.y = element_text(size = 14, color="black"),# Increase font size of pathway names
        axis.ticks = element_line(color = "black"),  # Add tick marks
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick mark length
        axis.line = element_line(color = "black"))
dev.off()