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

#Retrieving genes with significant temporal changes (from linear regression analysis)
CIS_Cor_Udegs <- read.xlsx("Supplementary_Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Cortex_Up")
CIS_Cor_Ddegs <- read.xlsx("Supplementary_Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Cortex_Down")

CIS_Int_Udegs <- read.xlsx("Supplementary_Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Interface_Up")
CIS_Int_Ddegs <- read.xlsx("Supplementary_Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Interface_Down")

CIS_Med_Udegs <- read.xlsx("Supplementary_Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Medulla_Up")
CIS_Med_Ddegs <- read.xlsx("Supplementary_Tables/Limma_cpmCIS_DEGs.xlsx", sheet = "CIS_Medulla_Down")

##### PLOTTING DISTRIBUTION OF LINEAR FEGRESSION SLOPES ########################
CIS_UAll<-rbind(CIS_Cor_Udegs, CIS_Int_Udegs, CIS_Med_Udegs)

CIS_DAll<-rbind(CIS_Cor_Ddegs, CIS_Int_Ddegs, CIS_Med_Ddegs)
#CIS_DAll$Slope<-abs(CIS_DAll$Slope)

CIS_AllSlope<-rbind(CIS_UAll,CIS_DAll)

  head(CIS_AllSlope)

data <- data.frame(
  value = c(log10(CIS_Cor_Udegs$Slope+1)/max(log10(abs(CIS_AllSlope$Slope)+1)),
            log10(CIS_Int_Udegs$Slope+1)/max(log10(abs(CIS_AllSlope$Slope)+1)), 
            log10(CIS_Med_Udegs$Slope+1)/max(log10(abs(CIS_AllSlope$Slope)+1)),
            log10(abs(CIS_Cor_Ddegs$Slope)+1)/max(log10(abs(CIS_AllSlope$Slope)+1)),
            log10(abs(CIS_Int_Ddegs$Slope)+1)/max(log10(abs(CIS_AllSlope$Slope)+1)), 
            log10(abs(CIS_Med_Ddegs$Slope)+1)/max(log10(abs(CIS_AllSlope$Slope)+1)) ),
  
  group = factor(c(rep("CIS_UCor", nrow(CIS_Cor_Udegs)), 
                   rep("CIS_UInt",  nrow(CIS_Int_Udegs)),
                   rep("CIS_UMed",  nrow(CIS_Med_Udegs)),
                   rep("CIS_DCor",  nrow(CIS_Cor_Ddegs)), 
                   rep("CIS_DInt",  nrow(CIS_Int_Ddegs)),
                   rep("CIS_DMed",  nrow(CIS_Med_Ddegs))  ))
)


data$group<-factor(data$group, levels=c("CIS_UCor", "CIS_UInt",  "CIS_UMed",
                                        "CIS_DCor", "CIS_DInt",  "CIS_DMed"))

# Plot the violin plot
#pdf("Figures/Figure2/pdfs/CIS_Slopes_Distribution.pdf", height=4, width=15)
ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7, width=1) + 
  geom_boxplot(width=0.05) +
  # Violin plot with transparency
  #geom_boxplot(width = 0.1, outlier.shape = TRUE, alpha = 0.9, color = "black") +    # Boxplot with black outline
  
  theme_minimal() +          # Use a clean theme
  labs(x = "Group", y = "Norm. log (abs((Slope)", title= "CIS Slope Distribution") +
  scale_fill_brewer(palette = "Set3") + 
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.ticks = element_line(color = "black"), 
        axis.ticks.length = unit(0.1, "cm")) +
  theme(legend.position = "none")# Optional color palette
 # dev.off()

