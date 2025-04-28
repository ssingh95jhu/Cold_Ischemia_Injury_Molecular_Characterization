#Clean environment
rm(list=ls(all.names=TRUE))
gc()

setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI/Manuscript")

########################## ADDING LIBRARY ######################################
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)

################################################################################
#### Extracting data from the stored tables:
## DEGS DETECTED IN DIFFERENT COMPARTMENTS OF THE KIDNEY
aki_cortex_degs <- read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Cortex")
rownames(aki_cortex_degs)<-aki_cortex_degs$Genes
head(aki_cortex_degs)

aki_interface_degs <- read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Interface")
rownames(aki_interface_degs)<-aki_interface_degs$Genes
head(aki_interface_degs)

aki_medulla_degs <- read.xlsx("Supplementary_Tables/Limma_cpmAKI2d_DEGs.xlsx", sheet = "AKI_Medulla")
rownames(aki_medulla_degs)<-aki_medulla_degs$Genes
head(aki_medulla_degs)


######## VISUALIZING COMMON GENES BETWEEN CIS COMPARTMENTS #####################
wb1<-createWorkbook()

####### 1.AKI CORTEX GENES vs AKI INTERFACE GENES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_genes<-intersect(rownames(aki_cortex_degs), rownames(aki_interface_degs) )
df<-data.frame(AKI_CORTEX_DEGS=aki_cortex_degs[comm_genes,]$Slope, AKI_INTERFACE_DEGS=aki_interface_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)

df$quadrant <- "green"
df[which(df$AKI_CORTEX_DEGS*df$AKI_INTERFACE_DEGS<0),]$quadrant<-"darkviolet"

#1.a. Pearson Correlation
Cor_Int_Corr<-cor.test(x=df$AKI_CORTEX_DEGS, y=df$AKI_INTERFACE_DEGS, method="pearson" )
Cor_Int_Corr <- data.frame(
  Statistic = Cor_Int_Corr$statistic,
  Correlation = Cor_Int_Corr$estimate,
  p_value = Cor_Int_Corr$p.value,
  Conf_Lower = Cor_Int_Corr$conf.int[1],
  Conf_Upper = Cor_Int_Corr$conf.int[2]
)

df_covary<-df[which(df$AKI_CORTEX_DEGS*df$AKI_INTERFACE_DEGS>0),]
df_covary<-data.frame(Gene=df_covary$gene, Cortex_Slope=df_covary$AKI_CORTEX_DEGS, Interface_Slope=df_covary$AKI_INTERFACE_DEGS)
head(df_covary)

df_disjoint<-df[which(df$AKI_CORTEX_DEGS*df$AKI_INTERFACE_DEGS<0),]
df_disjoint<-data.frame(Gene=df_disjoint$gene, Cortex_Slope=df_disjoint$AKI_CORTEX_DEGS, Interface_Slope=df_disjoint$AKI_INTERFACE_DEGS)
head(df_disjoint)

addWorksheet(wb1, "Cov_Cor_Int")
writeData(wb1, "Cov_Cor_Int", df_covary)

addWorksheet(wb1, "Dis_Cor_Int")
writeData(wb1, "Dis_Cor_Int", df_disjoint)

addWorksheet(wb1, "Corr_Cor_Int")
writeData(wb1, "Corr_Cor_Int", Cor_Int_Corr)

highlight_list<-c("Havcr1", "Lcn2", "Igfbp7","Spp1")

model <- lm(AKI_INTERFACE_DEGS ~ AKI_CORTEX_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

AKI_Cor_Int_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(AKI_Cor_Int_Fit)

addWorksheet(wb1, "R2_Cor_Int")
writeData(wb1, "R2_Cor_Int", AKI_Cor_Int_Fit)

slope <- coef(model)[2]  # Extract slope
intercept <- coef(model)[1]  # Extract intercept

pdf("Figures/Figure4S/pdfs/AKI_DEGs_Cortex_vs_Interface_32000x32000.pdf", height=5, width=7.5)
ggplot(df, aes(x = AKI_CORTEX_DEGS, y = AKI_INTERFACE_DEGS)) +
  # geom_point(alpha = 0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin = -2500, xmax = 2500, ymin = -2500, ymax = 2500), 
            fill = NA, color = 'red', linetype = 'solid', size = 0.5) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(AKI_CORTEX_DEGS) > 500 | abs(AKI_INTERFACE_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 30, force = 10, box.padding = 0.5) +  
  scale_x_continuous(limits = c(-21000, 11000)) +
  scale_y_continuous(limits = c(-21000, 11000)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()
  
#Inset plot for the above plot
pdf("Figures/Figure4S/pdfs/AKI_DEGs_Cortex_vs_Interface_5000x5000.pdf", height=5, width=7.5)
ggplot(df, aes(x = AKI_CORTEX_DEGS, y = AKI_INTERFACE_DEGS)) +
  # geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
  #data = \(d) d %>% filter((AKI_CORTEX_DEGS) <1000 | (AKI_INTERFACE_DEGS) <1000),
  data = \(d) d %>% filter((AKI_CORTEX_DEGS) >350 | (AKI_INTERFACE_DEGS) >350),
  #mapping = aes(label = gene, ), 
  mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
  max.overlaps = 45, box.padding=0.5) + 
  scale_x_continuous(limits = c(-2500, 2500)) +
  scale_y_continuous(limits = c(-2500, 2500))  + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1))
dev.off()

####### AKI CORTEX GENES vs AKI MEDULLA GENES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_genes<-intersect(rownames(aki_cortex_degs), rownames(aki_medulla_degs) )
df<-data.frame(AKI_CORTEX_DEGS=aki_cortex_degs[comm_genes,]$Slope, AKI_MEDULLA_DEGS=aki_medulla_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)

df$quadrant <- "green"
df[which(df$AKI_CORTEX_DEGS*df$AKI_MEDULLA_DEGS<0),]$quadrant<-"darkviolet"

Cor_Med_Corr<-cor.test(x=df$AKI_CORTEX_DEGS, y=df$AKI_MEDULLA_DEGS, method="pearson" )
Cor_Med_Corr<- data.frame(
  Statistic = Cor_Med_Corr$statistic,
  Correlation = Cor_Med_Corr$estimate,
  p_value = Cor_Med_Corr$p.value,
  Conf_Lower = Cor_Med_Corr$conf.int[1],
  Conf_Upper = Cor_Med_Corr$conf.int[2]
)

df_covary<-df[which(df$AKI_CORTEX_DEGS*df$AKI_MEDULLA_DEGS>0),]
df_covary<-data.frame(Gene=df_covary$gene, Cortex_Slope=df_covary$AKI_CORTEX_DEGS, Medulla_Slope=df_covary$AKI_MEDULLA_DEGS)
head(df_covary)

df_disjoint<-df[which(df$AKI_CORTEX_DEGS*df$AKI_MEDULLA_DEGS<0),]
df_disjoint<-data.frame(Gene=df_disjoint$gene, Cortex_Slope=df_disjoint$AKI_CORTEX_DEGS, Medulla_Slope=df_disjoint$AKI_MEDULLA_DEGS)
head(df_disjoint)

addWorksheet(wb1, "Cov_Cor_Med")
writeData(wb1, "Cov_Cor_Med", df_covary)

addWorksheet(wb1, "Dis_Cor_Med")
writeData(wb1, "Dis_Cor_Med", df_disjoint)

addWorksheet(wb1, "Corr_Cor_Med")
writeData(wb1, "Corr_Cor_Med", Cor_Med_Corr)

highlight_list<-c("Havcr1", "Lcn2", "Igfbp7","Spp1")

model <- lm(AKI_MEDULLA_DEGS ~ AKI_CORTEX_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

AKI_Cor_Med_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(AKI_Cor_Med_Fit)

addWorksheet(wb1, "R2_Cor_Med")
writeData(wb1, "R2_Cor_Med", AKI_Cor_Med_Fit)

slope <- coef(model)[2]  # Extract slope
intercept <- coef(model)[1]  # Extract intercept

pdf("Figures/Figure4S/pdfs/AKI_DEGs_Cortex_vs_Medulla_32000x32000.pdf", height=5, width=7.5)
ggplot(df, aes(x = AKI_CORTEX_DEGS, y = AKI_MEDULLA_DEGS)) +
  # geom_point(alpha = 0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin = -2500, xmax = 2500, ymin = -2500, ymax = 2500), 
            fill = NA, color = 'red', linetype = 'solid', size = 0.5) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(AKI_CORTEX_DEGS) > 500 | abs(AKI_MEDULLA_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 30, force = 10, box.padding = 0.5) +  
  scale_x_continuous(limits = c(-21000, 11000)) +
  scale_y_continuous(limits = c(-21000, 11000)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

#Inset plot for the above plot
pdf("Figures/Figure4S/pdfs/AKI_DEGs_Cortex_vs_Medulla_5000x5000.pdf", height=5, width=7.5)
ggplot(df, aes(x = AKI_CORTEX_DEGS, y = AKI_MEDULLA_DEGS)) +
  # geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter((AKI_CORTEX_DEGS) >300 | (AKI_MEDULLA_DEGS) >300),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 200, box.padding=0.5) + 
  scale_x_continuous(limits = c(-2500, 2500)) +
  scale_y_continuous(limits = c(-2500, 2500)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1))
dev.off()

####### AKI INTERFACE GENES vs AKI MEDULLA GENES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comm_genes<-intersect(rownames(aki_interface_degs), rownames(aki_medulla_degs) )
df<-data.frame(AKI_INTERFACE_DEGS=aki_interface_degs[comm_genes,]$Slope, AKI_MEDULLA_DEGS=aki_medulla_degs[comm_genes,]$Slope)
df$gene <- comm_genes
head(df)

df$quadrant <- "green"
df[which(df$AKI_INTERFACE_DEGS*df$AKI_MEDULLA_DEGS<0),]$quadrant<-"darkviolet"

Int_Med_Corr<-cor.test(x=df$AKI_INTERFACE_DEGS, y=df$AKI_MEDULLA_DEGS, method="pearson" )
Int_Med_Corr<- data.frame(
  Statistic = Int_Med_Corr$statistic,
  Correlation = Int_Med_Corr$estimate,
  p_value = Int_Med_Corr$p.value,
  Conf_Lower = Int_Med_Corr$conf.int[1],
  Conf_Upper = Int_Med_Corr$conf.int[2]
)

df_covary<-df[which(df$AKI_INTERFACE_DEGS*df$AKI_MEDULLA_DEGS>0),]
df_covary<-data.frame(Gene=df_covary$gene, Interface_Slope=df_covary$AKI_INTERFACE_DEGS, Medulla_Slope=df_covary$AKI_MEDULLA_DEGS)
head(df_covary)

df_disjoint<-df[which(df$AKI_INTERFACE_DEGS*df$AKI_MEDULLA_DEGS<0),]
df_disjoint<-data.frame(Gene=df_disjoint$gene, Interface_Slope=df_disjoint$AKI_INTERFACE_DEGS, Medulla_Slope=df_disjoint$AKI_MEDULLA_DEGS)
head(df_disjoint)

addWorksheet(wb1, "Cov_Int_Med")
writeData(wb1, "Cov_Int_Med", df_covary)

addWorksheet(wb1, "Dis_Int_Med")
writeData(wb1, "Dis_Int_Med", df_disjoint)

addWorksheet(wb1, "Corr_Int_Med")
writeData(wb1, "Corr_Int_Med", Int_Med_Corr)


highlight_list<-c("Havcr1", "Lcn2", "Igfbp7","Spp1")

model <- lm(AKI_MEDULLA_DEGS ~ AKI_INTERFACE_DEGS, data = df)
summary<-summary(model)
R2<-summary$r.squared # Look for R-squared
p<-summary$fstatistic
pvalue<-pf(p[1], p[2], p[3], lower.tail = FALSE)
print(pvalue)

AKI_Int_Med_Fit<-data.frame(R2=R2, pvalue=pvalue)
print(AKI_Int_Med_Fit)

addWorksheet(wb1, "R2_Int_Med")
writeData(wb1, "R2_Int_Med", AKI_Int_Med_Fit)

slope <- coef(model)[2]  # Extract slope
intercept <- coef(model)[1]  # Extract intercept

pdf("Figures/Figure4S/pdfs/AKI_DEGs_Interface_vs_Medulla_32000x32000.pdf", height=5, width=7.5)
ggplot(df, aes(x = AKI_INTERFACE_DEGS, y = AKI_MEDULLA_DEGS)) +
  #geom_point(alpha = 0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_rect(aes(xmin = -2500, xmax = 2500, ymin = -2500, ymax = 2500), 
            fill = NA, color = 'red', linetype = 'solid', size = 0.5) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter(abs(AKI_INTERFACE_DEGS) > 500 | abs(AKI_MEDULLA_DEGS) > 500),
    mapping = aes(label = gene, color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 30, force = 10, box.padding = 0.5
  ) +  
  scale_x_continuous(limits = c(-21000, 11000)) +
  scale_y_continuous(limits = c(-21000, 11000)) + 
  scale_color_identity() + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

pdf("Figures/Figure4S/pdfs/AKI_DEGs_Interface_vs_Medulla_5000x5000.pdf", height=5, width=7.5)
ggplot(df, aes(x = AKI_INTERFACE_DEGS, y = AKI_MEDULLA_DEGS)) +
  #geom_point(alpha=0.2) +
  geom_point(aes(color=quadrant),alpha=0.8, shape=16) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 0.5) +  # Pearson correlation line
  geom_text_repel(
    data = \(d) d %>% filter((AKI_INTERFACE_DEGS) >300 | (AKI_MEDULLA_DEGS) >300),
    mapping = aes(label = gene,  color = ifelse(gene %in% highlight_list, "red", "black")), 
    max.overlaps = 200, box.padding=0.5) + 
  scale_x_continuous(limits = c(-2500, 2500)) +
  scale_y_continuous(limits = c(-2500, 2500)) + 
  scale_color_identity() + theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
        panel.border = element_rect(color = "red", fill = NA, size = 1))
dev.off()

saveWorkbook(wb1,"Supplementary_Tables/Fig4S_AKI_Compartmental_Covary_Disjoint_DEGs.xlsx", overwrite=TRUE)

