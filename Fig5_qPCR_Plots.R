#Clean environment
rm(list=ls(all.names=TRUE))
gc()

library(ggplot2)
library(ggsci)

workdir <- setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI")

#### Pink1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g <- 'Pink1'
data <- read.csv(paste0(workdir, "/Ryo_Data/", g, ".csv"))
#data$Time<-gsub("h","",data$Time)
data$Region <- factor(data$Region, levels=c("Cortex", "Medulla"))
data$Time <- factor(data$Time, levels=c("0h", "12h","24h","48h"))
fig_dir<-paste0(workdir, "/Manuscript/")

#Plotted on the same plot
pdf(paste0(fig_dir,"/Figures/Figure5/pdfs/qPCR_Pink1.pdf"), height=4, width=6)
ggplot(data, aes(x = Time, y = RQ, color = Region, group = Region, fill = Region)) +  
  geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
  geom_smooth(method = 'lm', se = TRUE) +  
  scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
  scale_color_aaas(name = "") +  
  scale_y_continuous(limits=c(-1,6)) +
  labs(x= "Cold Ischemia Time", y = "Relative mRNA Fold Change", title = g) +  
  theme_minimal() +  
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(size = 12, face = "bold.italic", hjust = 0.5),  
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),  
        axis.text.x = element_text(size = 10, color ="black"), axis.text.y = element_text(size = 10, color = "black"),  
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"), 
        legend.text = element_text(size = 10))
dev.off()

