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
library(ggrepel)
library(dplyr)
library(patchwork)
library(cowplot)
library(viridis) #color blind friendly
library(RColorBrewer) #for RColorbrewer

################################################################################
#################### NORMALIZE CTRL DATASET ####################################

## limma performs better with gaussian distributed
ctrl1$mat <- MERINGUE::normalizeCounts(ctrl1$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
ctrl2$mat <- MERINGUE::normalizeCounts(ctrl2$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
ctrl3$mat <- MERINGUE::normalizeCounts(ctrl3$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
ctrl4$mat <- MERINGUE::normalizeCounts(ctrl4$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)

## easier for plotting
ctrl1$mat_notlog <- MERINGUE::normalizeCounts(ctrl1$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl2$mat_notlog <- MERINGUE::normalizeCounts(ctrl2$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl3$mat_notlog <- MERINGUE::normalizeCounts(ctrl3$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
ctrl4$mat_notlog <- MERINGUE::normalizeCounts(ctrl4$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

#################### NORMALIZE IRL DATASET ####################################

## limma performs better with gaussian distributed
irl1$mat <- MERINGUE::normalizeCounts(irl1$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
irl2$mat <- MERINGUE::normalizeCounts(irl2$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
irl3$mat <- MERINGUE::normalizeCounts(irl3$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
irl4$mat <- MERINGUE::normalizeCounts(irl4$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)

## easier for plotting
irl1$mat_notlog <- MERINGUE::normalizeCounts(irl1$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl2$mat_notlog <- MERINGUE::normalizeCounts(irl2$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl3$mat_notlog <- MERINGUE::normalizeCounts(irl3$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
irl4$mat_notlog <- MERINGUE::normalizeCounts(irl4$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

############################ normalize #########################################

## limma performs better with gaussian distributed
CIS_0h$mat <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
CIS_12h$mat <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
CIS_24h$mat <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
CIS_48h$mat <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)

## easier for plotting
CIS_0h$mat_notlog <- MERINGUE::normalizeCounts(CIS_0h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_12h$mat_notlog <- MERINGUE::normalizeCounts(CIS_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_24h$mat_notlog <- MERINGUE::normalizeCounts(CIS_24h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
CIS_48h$mat_notlog <- MERINGUE::normalizeCounts(CIS_48h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

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

################# normalize AKI DATASET (CPM) ##################################

AKI_sham$mat <- MERINGUE::normalizeCounts(AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
AKI_4h$mat <- MERINGUE::normalizeCounts(AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
AKI_12h$mat <- MERINGUE::normalizeCounts(AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
AKI_2d$mat <- MERINGUE::normalizeCounts(AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)
AKI_6w$mat <- MERINGUE::normalizeCounts(AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),], log=TRUE)

AKI_sham$mat_notlog <- MERINGUE::normalizeCounts(AKI_sham$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_4h$mat_notlog <- MERINGUE::normalizeCounts(AKI_4h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_12h$mat_notlog <- MERINGUE::normalizeCounts(AKI_12h$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_2d$mat_notlog <- MERINGUE::normalizeCounts(AKI_2d$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)
AKI_6w$mat_notlog <- MERINGUE::normalizeCounts(AKI_6w$gexp[unique(rownames(CIS_0h$gexp)),], log=FALSE)

##############################################################################
#Extracting the compartment specific spots for CIS and AKI:

cortex <- readRDS('AKI-CIS-irl-ctrl_Cortex_spots.rds')
interface <- readRDS('AKI-CIS-irl-ctrl_Interface_spots.rds')
medulla <- readRDS('AKI-CIS-irl-ctrl_Medulla_spots.rds')
#all_spots<-readRDS('AKI-CIS-irl-ctrl_All_spots.rds')

#1.All_Cortex_Spots (CTRL)
CTRL1_cortex<-names(cortex[grep("^CTRL1", names(cortex))])
CTRL2_cortex<-names(cortex[grep("^CTRL2", names(cortex))])
CTRL3_cortex<-names(cortex[grep("^CTRL3", names(cortex))])
CTRL4_cortex<-names(cortex[grep("^CTRL4", names(cortex))])

#2.All_Interface_Spots (CTRL)
CTRL1_interface<-names(interface[grep("^CTRL1", names(interface))])
CTRL2_interface<-names(interface[grep("^CTRL2", names(interface))])
CTRL3_interface<-names(interface[grep("^CTRL3", names(interface))])
CTRL4_interface<-names(interface[grep("^CTRL4", names(interface))])

#3.All_Medulla_Spots (CTRL)
CTRL1_medulla<-names(medulla[grep("^CTRL1", names(medulla))])
CTRL2_medulla<-names(medulla[grep("^CTRL2", names(medulla))])
CTRL3_medulla<-names(medulla[grep("^CTRL3", names(medulla))])
CTRL4_medulla<-names(medulla[grep("^CTRL4", names(medulla))])

#1.All_Cortex_Spots (IRL)
IRL1_cortex<-names(cortex[grep("^IR1", names(cortex))])
IRL2_cortex<-names(cortex[grep("^IR2", names(cortex))])
IRL3_cortex<-names(cortex[grep("^IR3", names(cortex))])
IRL4_cortex<-names(cortex[grep("^IR4", names(cortex))])

#2.All_Interface_Spots (IRL)
IRL1_interface<-names(interface[grep("^IR1", names(interface))])
IRL2_interface<-names(interface[grep("^IR2", names(interface))])
IRL3_interface<-names(interface[grep("^IR3", names(interface))])
IRL4_interface<-names(interface[grep("^IR4", names(interface))])

#3.All_Medulla_Spots (IRL)
IRL1_medulla<-names(medulla[grep("^IR1", names(medulla))])
IRL2_medulla<-names(medulla[grep("^IR2", names(medulla))])
IRL3_medulla<-names(medulla[grep("^IR3", names(medulla))])
IRL4_medulla<-names(medulla[grep("^IR4", names(medulla))])


#1.All_Cortex_Spots (CIS)
CIS_cortex.0h<-names(cortex[grep("^CIS_0h", names(cortex))])
CIS_cortex.12h<-names(cortex[grep("^CIS_12h", names(cortex))])
CIS_cortex.24h<-names(cortex[grep("^CIS_24h", names(cortex))])
CIS_cortex.48h<-names(cortex[grep("^CIS_48h", names(cortex))])

#2.All_Interface_Spots (CIS)
CIS_interface.0h<-names(interface[grep("^CIS_0h", names(interface))])
CIS_interface.12h<-names(interface[grep("^CIS_12h", names(interface))])
CIS_interface.24h<-names(interface[grep("^CIS_24h", names(interface))])
CIS_interface.48h<-names(interface[grep("^CIS_48h", names(interface))])

#3.All_Medulla_Spots (CIS)
CIS_medulla.0h<-names(medulla[grep("^CIS_0h", names(medulla))])
CIS_medulla.12h<-names(medulla[grep("^CIS_12h", names(medulla))])
CIS_medulla.24h<-names(medulla[grep("^CIS_24h", names(medulla))])
CIS_medulla.48h<-names(medulla[grep("^CIS_48h", names(medulla))])


#1.All_Cortex_Spots (AKI)
AKI_cortex.sham<-names(cortex[grep("^AKI_sham", names(cortex))])
AKI_cortex.4h<-names(cortex[grep("^AKI_4h", names(cortex))])
AKI_cortex.12h<-names(cortex[grep("^AKI_12h", names(cortex))])
AKI_cortex.2d<-names(cortex[grep("^AKI_2d", names(cortex))])
AKI_cortex.6w<-names(cortex[grep("^AKI_6w", names(cortex))])

#2.All_Interface_Spots (AKI)
AKI_interface.sham<-names(interface[grep("^AKI_sham", names(interface))])
AKI_interface.4h<-names(interface[grep("^AKI_4h", names(interface))])
AKI_interface.12h<-names(interface[grep("^AKI_12h", names(interface))])
AKI_interface.2d<-names(interface[grep("^AKI_2d", names(interface))])
AKI_interface.6w<-names(interface[grep("^AKI_6w", names(interface))])

#3.All_Medulla_Spots (AKI)
AKI_medulla.sham<-names(medulla[grep("^AKI_sham", names(medulla))])
AKI_medulla.4h<-names(medulla[grep("^AKI_4h", names(medulla))])
AKI_medulla.12h<-names(medulla[grep("^AKI_12h", names(medulla))])
AKI_medulla.2d<-names(medulla[grep("^AKI_2d", names(medulla))])
AKI_medulla.6w<-names(medulla[grep("^AKI_6w", names(medulla))])

####################### PLOTTING AND VISUALIZATION #############################
##################### PLOTTING (CIS ONLY)
g<-'Tfeb'#'Rack1'#'Ccn2' #Angpt1'#'Cntf'
#'Clcnka'#'Aqp1'#'Cat'#'Hif1a'#'Ppara'#Egf'#'Esrra'#'Cpt1a'#'Cd36'#'Il2rg'#Gsr'#'Irf1'#'Sod1'#'Cat'#'Sod2'#'Atp5a1'#'Sod3'

overall_min <- min(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                   CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,])

overall_max <- quantile(c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                          CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), 0.95)

color_scale <- scale_color_gradientn(colors=c("darkgray","yellow","red"),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1), 
                                     limits = c(overall_min, overall_max), oob = scales::squish)

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

df <- data.frame(CIS_0h$pos, gexp=CIS_0h$mat_notlog[g,])
g1 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme  + ggtitle('CIS 0h') 
df <- data.frame(CIS_12h$pos, gexp=CIS_12h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1)  + color_scale  + my_theme + ggtitle('CIS 12h')
df <- data.frame(CIS_24h$pos, gexp=CIS_24h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme + ggtitle('CIS 24h')
df <- data.frame(CIS_48h$pos, gexp=CIS_48h$mat_notlog[g,])
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('CIS 48h')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + my_theme + ggtitle('CIS 48h')

comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))

df1 <- data.frame(t = c(0, 12, 24, 48),
                  gexp = c(mean(CIS_0h$mat_notlog[g,CIS_cortex.0h]), mean(CIS_12h$mat_notlog[g,CIS_cortex.12h]), 
                           mean(CIS_24h$mat_notlog[g,CIS_cortex.24h]), mean(CIS_48h$mat_notlog[g,CIS_cortex.48h])), group = 'Cortex')
df2 <- data.frame(t = c(0, 12, 24, 48),
                  gexp = c(mean(CIS_0h$mat_notlog[g,CIS_interface.0h]), mean(CIS_12h$mat_notlog[g,CIS_interface.12h]), 
                           mean(CIS_24h$mat_notlog[g,CIS_interface.24h]), mean(CIS_48h$mat_notlog[g,CIS_interface.48h])), group = 'Interface')

df3 <- data.frame(t = c(0, 12, 24, 48),
                  gexp = c(mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                           mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h])), group = 'Medulla')

df4 <- data.frame(t = c(0, 12, 24, 48),
                  gexp = c(mean(CIS_0h$mat_notlog[g,]), mean(CIS_12h$mat_notlog[g,]), 
                           mean(CIS_24h$mat_notlog[g,]), mean(CIS_48h$mat_notlog[g,])), group = 'Global')

overall_min0<-min(mean(CIS_0h$mat_notlog[g,]), mean(CIS_12h$mat_notlog[g,]), 
                  mean(CIS_24h$mat_notlog[g,]), mean(CIS_48h$mat_notlog[g,]),
                  mean(CIS_0h$mat_notlog[g,CIS_cortex.0h]), mean(CIS_12h$mat_notlog[g,CIS_cortex.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_cortex.24h]), mean(CIS_48h$mat_notlog[g,CIS_cortex.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_interface.0h]), mean(CIS_12h$mat_notlog[g,CIS_interface.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_interface.24h]), mean(CIS_48h$mat_notlog[g,CIS_interface.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h]) )

overall_max0<-max(mean(CIS_0h$mat_notlog[g,]), mean(CIS_12h$mat_notlog[g,]), 
                  mean(CIS_24h$mat_notlog[g,]), mean(CIS_48h$mat_notlog[g,]),
                  mean(CIS_0h$mat_notlog[g,CIS_cortex.0h]), mean(CIS_12h$mat_notlog[g,CIS_cortex.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_cortex.24h]), mean(CIS_48h$mat_notlog[g,CIS_cortex.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_interface.0h]), mean(CIS_12h$mat_notlog[g,CIS_interface.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_interface.24h]), mean(CIS_48h$mat_notlog[g,CIS_interface.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h]) )

df<-rbind(df1,df2,df3, df4)
head(df)
df$group<-factor(df$group, levels=c("Global","Cortex","Interface","Medulla"))
# gline <- ggplot(df, aes(x = t, y = gexp, col=group)) + geom_line() + geom_point() + 
#          ylim(overall_min0, overall_max0) + theme_light() +
#          theme(axis.text.x = element_text(face = "bold", color="black"),  
#                axis.text.y = element_text(face = "bold", color="black"),
#                panel.border = element_rect(color = "black", fill = NA, size = 1),
#                legend.text = element_text(face = "bold"),
#                panel.grid.minor = element_blank() ) +
#          xlab("time (hours)") + ylab("Avg. gene exp. (counts)")


gline <- ggplot(df, aes(x = t, y = gexp, col=group)) + geom_line() + geom_point() + 
         scale_color_manual(values=c("Global"="black","Cortex"="red","Interface"="blue","Medulla"="green")) +
         ylim(overall_min0, overall_max0) + theme_light() +
         theme(axis.text.x = element_text(color="black", size=11),  
         axis.text.y = element_text(color="black",size=11),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=12),
         panel.border = element_rect(color = "black", fill = NA, size = 1),
         legend.text = element_text(color="black", size=12),
         panel.grid.minor = element_blank(),
         axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.1, "cm")) +
         xlab("time (hours)") + ylab("Avg. gene exp. (counts)")

library(gridExtra)
#grid.arrange(g1,g2,g3,g4,gline, ncol=5, top=g)
pdf("Figures/FigureY/pdfs/CIS_Tfeb_AllTimepoints.pdf", height=3, width=15)
grid.arrange(comb_plot,gline, ncol=2, widths=c(16, 5), top=g)
dev.off()

###############################################################################
##################### PLOTTING (AKI ONLY)
g<-'Lcn2'#'Egf'#'Cpt1a'#'Cd36'#'Il2rg'#'Gsr'#'Irf1'#'B2m'#'Sod1'#'Cat'#'Dlst'#Sod2'#'Atp5a1'#'Sod3'

overall_min <- min(AKI_sham$mat_notlog[g,], AKI_4h$mat_notlog[g,], 
                   AKI_12h$mat_notlog[g,], AKI_2d$mat_notlog[g,])

overall_max <- quantile(c(AKI_sham$mat_notlog[g,], AKI_4h$mat_notlog[g,], 
                          AKI_12h$mat_notlog[g,], AKI_2d$mat_notlog[g,]), 0.95)

color_scale <- scale_color_gradientn(colors=c("darkgray","yellow","red"),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1), 
                                     limits = c(overall_min, overall_max), oob = scales::squish)

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

df <- data.frame(AKI_sham$pos, gexp=AKI_sham$mat_notlog[g,])
g1 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('AKI sham')
df <- data.frame(AKI_4h$pos, gexp=AKI_4h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme + ggtitle('AKI 4h')
df <- data.frame(AKI_12h$pos, gexp=AKI_12h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('AKI 12h')
df <- data.frame(AKI_2d$pos, gexp=AKI_2d$mat_notlog[g,])
g4 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale +  theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('AKI 2d')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale  + my_theme + ggtitle('AKI 2d')

comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))


df1 <- data.frame(t = c(0, 4, 12, 48),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d])), group = 'Cortex')

df2 <- data.frame(t = c(0, 4, 12, 48),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d])), group = 'Interface')

df3 <- data.frame(t = c(0, 4, 12, 48),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d])), group = 'Medulla')

df4 <- data.frame(t = c(0, 4, 12, 48),
                  gexp = c(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                           mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,])), group = 'Global')

overall_min0<-min(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                  mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,]),
                  mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d]),
                  mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d]),
                  mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d]) )

overall_max0<-max(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                  mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,]),
                  mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d]),
                  mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d]),
                  mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d]) )

df<-rbind(df1,df2,df3, df4)
head(df)
df$group<-factor(df$group, levels=c("Global","Cortex","Interface","Medulla"))

gline <- ggplot(df, aes(x = t, y = gexp, col=group)) + geom_line() + geom_point() + 
  scale_color_manual(values=c("Global"="black","Cortex"="red","Interface"="blue","Medulla"="green")) +
  ylim(overall_min0, overall_max0) + theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(0.1, "cm")) +
  xlab("time (hours)") + ylab("Avg. gene exp. (counts)")

library(gridExtra)
#grid.arrange(g1,g2,g3,g4,gline, ncol=5, top=g)
pdf("Figures/FigureX/pdfs/AKI_Lcn2_AllTimepoints.pdf", height=3, width=15)
grid.arrange(comb_plot,gline, ncol=2, widths=c(16, 5), top=g)
dev.off()

############## AKI ONLY (Sham - 6 weeks)

g<-'Cox6a1'#'Egf'#'Cpt1a'#'Cd36'#'Il2rg'#'Gsr'#'Irf1'#'B2m'#'Sod1'#'Cat'#'Dlst'#Sod2'#'Atp5a1'#'Sod3'

overall_min <- min(AKI_sham$mat_notlog[g,], AKI_4h$mat_notlog[g,], 
                   AKI_12h$mat_notlog[g,], AKI_2d$mat_notlog[g,])

overall_max <- quantile(c(AKI_sham$mat_notlog[g,], AKI_4h$mat_notlog[g,], 
                          AKI_12h$mat_notlog[g,], AKI_2d$mat_notlog[g,]), 0.95)

color_scale <- scale_color_gradientn(colors=c("darkgray","yellow","red"),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1), 
                                     limits = c(overall_min, overall_max), oob = scales::squish)

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)



my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

df <- data.frame(AKI_sham$pos, gexp=AKI_sham$mat_notlog[g,])
g1 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('AKI sham')
df <- data.frame(AKI_4h$pos, gexp=AKI_4h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme + ggtitle('AKI 4h')
df <- data.frame(AKI_12h$pos, gexp=AKI_12h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('AKI 12h')
df <- data.frame(AKI_2d$pos, gexp=AKI_2d$mat_notlog[g,])
g4 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('AKI 2d')
df <- data.frame(AKI_6w$pos, gexp=AKI_6w$mat_notlog[g,])
g5 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale +  theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('AKI 6w')


# Extract the legend from one plot
legend <- get_legend(g5)
g5 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale  + my_theme + ggtitle('AKI 6w')

comb_plot<-plot_grid(g1, g2, g3, g4, g5, legend, ncol = 6, rel_widths = c(1, 1, 1, 1, 1, 0.25))


df1 <- data.frame(t = c(0, 4, 12, 48, 1008),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d]), 
                           mean(AKI_6w$mat_notlog[g,AKI_cortex.6w])), group = 'Cortex')

df2 <- data.frame(t = c(0, 4, 12, 48, 1008),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d]), 
                           mean(AKI_6w$mat_notlog[g,AKI_interface.6w])),group = 'Interface')

df3 <- data.frame(t = c(0, 4, 12, 48, 1008),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d]), 
                           mean(AKI_6w$mat_notlog[g,AKI_medulla.6w])),group = 'Medulla')

df4 <- data.frame(t = c(0, 4, 12, 48, 1008),
                  gexp = c(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                           mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,]), 
                           mean(AKI_6w$mat_notlog[g,])), group = 'Global')

overall_min0<-min(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                  mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,]), mean(AKI_6w$mat_notlog[g,]),
                  mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d]), mean(AKI_6w$mat_notlog[g,AKI_cortex.6w]),
                  mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d]), mean(AKI_6w$mat_notlog[g,AKI_interface.6w]),
                  mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d]), mean(AKI_6w$mat_notlog[g,AKI_medulla.6w]) )

overall_max0<-max(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                  mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,]), mean(AKI_6w$mat_notlog[g,]),
                  mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d]), mean(AKI_6w$mat_notlog[g,AKI_cortex.6w]),
                  mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d]), mean(AKI_6w$mat_notlog[g,AKI_interface.6w]),
                  mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d]), mean(AKI_6w$mat_notlog[g,AKI_medulla.6w]) )

df<-rbind(df1,df2,df3, df4)
head(df)
df$group<-factor(df$group, levels=c("Global","Cortex","Interface","Medulla"))

gline <- ggplot(df, aes(x = log10(t+1), y = gexp, col=group)) + geom_line() + geom_point() + 
  scale_color_manual(values=c("Global"="black","Cortex"="red","Interface"="blue","Medulla"="green")) +
  ylim(overall_min0, overall_max0) + theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank() ) +
  xlab("time (hours)") + ylab("Avg. gene exp. (counts)")

library(gridExtra)
#grid.arrange(g1,g2,g3,g4,gline, ncol=5, top=g)
#pdf("Figures/Figure3/pdfs/AKI_Aldoa_AllTimepoints.pdf", height=3, width=15)
grid.arrange(comb_plot,gline, ncol=2, widths=c(20, 5), top=g)
#dev.off()




##################### PLOTTING (CTRL)
                                   

###################### PLOTTING IRL ONLY #######################################
g ='Mdh1'#'Srebf2'#'Fgg'

overall_min <- min(irl1$mat_notlog[g,], irl2$mat_notlog[g,], irl3$mat_notlog[g,], irl4$mat_notlog[g,])

overall_max <- quantile(c(irl1$mat_notlog[g,], irl2$mat_notlog[g,], irl3$mat_notlog[g,], irl4$mat_notlog[g,]), 0.95)


color_scale <- scale_color_gradientn(colors=c("darkgray","yellow","red"),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

df <- data.frame(irl1$pos, gexp=irl1$mat_notlog[g,])
g1 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('IRL1')
df <- data.frame(irl2$pos, gexp=irl2$mat_notlog[g,])
g2 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme + ggtitle('IRL2')
df <- data.frame(irl3$pos, gexp=irl3$mat_notlog[g,])
g3 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('IRL3')
df <- data.frame(irl4$pos, gexp=irl4$mat_notlog[g,])
g4 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale +  theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('AKI 2d')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale  + my_theme + ggtitle('IRL4')

comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))

df1 <- data.frame(t = c(1, 2, 3, 4),
                  gexp = c(mean(irl1$mat_notlog[g,IRL1_cortex]), mean(irl2$mat_notlog[g,IRL2_cortex]), 
                           mean(irl3$mat_notlog[g,IRL3_cortex]), mean(irl4$mat_notlog[g,IRL4_cortex])), group = 'Cortex')

df2 <- data.frame(t = c(1, 2, 3, 4),
                  gexp = c(mean(irl1$mat_notlog[g,IRL1_interface]), mean(irl2$mat_notlog[g,IRL2_interface]), 
                           mean(irl3$mat_notlog[g,IRL3_interface]), mean(irl4$mat_notlog[g,IRL4_interface])), group = 'Interface')

df3 <- data.frame(t = c(1, 2, 3, 4),
                  gexp = c(mean(irl1$mat_notlog[g,IRL1_medulla]), mean(irl2$mat_notlog[g,IRL2_medulla]), 
                           mean(irl3$mat_notlog[g,IRL3_medulla]), mean(irl4$mat_notlog[g,IRL4_medulla])), group = 'Medulla')

df4 <- data.frame(t = c(1, 2, 3, 4),
                  gexp = c(mean(irl1$mat_notlog[g,]), mean(irl2$mat_notlog[g,]), 
                           mean(irl3$mat_notlog[g,]), mean(irl4$mat_notlog[g,])), group = 'Global')

overall_min0<-min(mean(irl1$mat_notlog[g,]), mean(irl2$mat_notlog[g,]), 
                  mean(irl3$mat_notlog[g,]), mean(irl4$mat_notlog[g,]),
                  mean(irl1$mat_notlog[g,IRL1_cortex]), mean(irl2$mat_notlog[g,IRL2_cortex]), 
                  mean(irl3$mat_notlog[g,IRL3_cortex]), mean(irl4$mat_notlog[g,IRL4_cortex]),
                  mean(irl1$mat_notlog[g,IRL1_interface]), mean(irl2$mat_notlog[g,IRL2_interface]), 
                  mean(irl3$mat_notlog[g,IRL3_interface]), mean(irl4$mat_notlog[g,IRL4_interface]),
                  mean(irl1$mat_notlog[g,IRL1_medulla]), mean(irl2$mat_notlog[g,IRL2_medulla]), 
                  mean(irl3$mat_notlog[g,IRL3_medulla]), mean(irl4$mat_notlog[g,IRL4_medulla]) )

overall_max0<-max(mean(irl1$mat_notlog[g,]), mean(irl2$mat_notlog[g,]), 
                  mean(irl3$mat_notlog[g,]), mean(irl4$mat_notlog[g,]),
                  mean(irl1$mat_notlog[g,IRL1_cortex]), mean(irl2$mat_notlog[g,IRL2_cortex]), 
                  mean(irl3$mat_notlog[g,IRL3_cortex]), mean(irl4$mat_notlog[g,IRL4_cortex]),
                  mean(irl1$mat_notlog[g,IRL1_interface]), mean(irl2$mat_notlog[g,IRL2_interface]), 
                  mean(irl3$mat_notlog[g,IRL3_interface]), mean(irl4$mat_notlog[g,IRL4_interface]),
                  mean(irl1$mat_notlog[g,IRL1_medulla]), mean(irl2$mat_notlog[g,IRL2_medulla]), 
                  mean(irl3$mat_notlog[g,IRL3_medulla]), mean(irl4$mat_notlog[g,IRL4_medulla]) )

df<-rbind(df1,df2,df3, df4)
head(df)
df$group<-factor(df$group, levels=c("Global","Cortex","Interface","Medulla"))

gline <- ggplot(df, aes(x = t, y = gexp, col=group)) + geom_point() + 
  scale_color_manual(values=c("Global"="black","Cortex"="red","Interface"="blue","Medulla"="green")) +
  ylim(overall_min0, overall_max0) + theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank() ) +
  xlab("Kidney_ID") + ylab("Avg. gene exp. (counts)")

library(gridExtra)
pdf("Figures/FigureX/pdfs/AKI24_Spp1_AllTimepoints.pdf", height=3, width=15)
grid.arrange(comb_plot,gline, ncol=2, widths=c(16, 5), top=g)
dev.off()

##################### PLOTTING (CTRL vs CIS vs AKI) ############################
g <-'Ppargc1a'#'Srebf2' #Hsp90b1'#'Hspa5'
overall_min <- min(ctrl1$mat_notlog[g,], ctrl2$mat_notlog[g,], ctrl3$mat_notlog[g,], ctrl4$mat_notlog[g,])

overall_max <- quantile(c(ctrl1$mat_notlog[g,], ctrl2$mat_notlog[g,], ctrl3$mat_notlog[g,], ctrl4$mat_notlog[g,]), 0.95)


color_scale <- scale_color_gradientn(colors=c("darkgray","yellow","red"),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)


my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

df <- data.frame(ctrl1$pos, gexp=ctrl1$mat_notlog[g,])
g1 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('CTRL1')
df <- data.frame(ctrl2$pos, gexp=ctrl2$mat_notlog[g,])
g2 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme + ggtitle('CTRL2')
df <- data.frame(ctrl3$pos, gexp=ctrl3$mat_notlog[g,])
g3 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('CTRL3')
df <- data.frame(ctrl4$pos, gexp=ctrl4$mat_notlog[g,])
g4 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale +  theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('AKI 2d')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=x, y=y, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale  + my_theme + ggtitle('CTRL4')

ctrl_comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))

df1 <- data.frame(t = c(1, 2, 3, 4),
                  gexp = c(mean(ctrl1$mat_notlog[g,CTRL1_cortex]), mean(ctrl2$mat_notlog[g,CTRL2_cortex]), 
                           mean(ctrl3$mat_notlog[g,CTRL3_cortex]), mean(ctrl4$mat_notlog[g,CTRL4_cortex])), group = 'Cortex')

df2 <- data.frame(t = c(1, 2, 3, 4),
                  gexp = c(mean(ctrl1$mat_notlog[g,CTRL1_interface]), mean(ctrl2$mat_notlog[g,CTRL2_interface]), 
                           mean(ctrl3$mat_notlog[g,CTRL3_interface]), mean(ctrl4$mat_notlog[g,CTRL4_interface])), group = 'Interface')

df3 <- data.frame(t = c(1, 2, 3, 4),
                  gexp = c(mean(ctrl1$mat_notlog[g,CTRL1_medulla]), mean(ctrl2$mat_notlog[g,CTRL2_medulla]), 
                           mean(ctrl3$mat_notlog[g,CTRL3_medulla]), mean(ctrl4$mat_notlog[g,CTRL4_medulla])), group = 'Medulla')

df4 <- data.frame(t = c(1, 2, 3, 4),
                  gexp = c(mean(ctrl1$mat_notlog[g,]), mean(ctrl2$mat_notlog[g,]), 
                           mean(ctrl3$mat_notlog[g,]), mean(ctrl4$mat_notlog[g,])), group = 'Global')

overall_min0<-min(mean(ctrl1$mat_notlog[g,]), mean(ctrl2$mat_notlog[g,]), 
                  mean(ctrl3$mat_notlog[g,]), mean(ctrl4$mat_notlog[g,]),
                  mean(ctrl1$mat_notlog[g,CTRL1_cortex]), mean(ctrl2$mat_notlog[g,CTRL2_cortex]), 
                  mean(ctrl3$mat_notlog[g,CTRL3_cortex]), mean(ctrl4$mat_notlog[g,CTRL4_cortex]),
                  mean(ctrl1$mat_notlog[g,CTRL1_interface]), mean(ctrl2$mat_notlog[g,CTRL2_interface]), 
                  mean(ctrl3$mat_notlog[g,CTRL3_interface]), mean(ctrl4$mat_notlog[g,CTRL4_interface]),
                  mean(ctrl1$mat_notlog[g,CTRL1_medulla]), mean(ctrl2$mat_notlog[g,CTRL2_medulla]), 
                  mean(ctrl3$mat_notlog[g,CTRL3_medulla]), mean(ctrl4$mat_notlog[g,CTRL4_medulla]) )

overall_max0<-max(mean(ctrl1$mat_notlog[g,]), mean(ctrl2$mat_notlog[g,]), 
                  mean(ctrl3$mat_notlog[g,]), mean(ctrl4$mat_notlog[g,]),
                  mean(ctrl1$mat_notlog[g,CTRL1_cortex]), mean(ctrl2$mat_notlog[g,CTRL2_cortex]), 
                  mean(ctrl3$mat_notlog[g,CTRL3_cortex]), mean(ctrl4$mat_notlog[g,CTRL4_cortex]),
                  mean(ctrl1$mat_notlog[g,CTRL1_interface]), mean(ctrl2$mat_notlog[g,CTRL2_interface]), 
                  mean(ctrl3$mat_notlog[g,CTRL3_interface]), mean(ctrl4$mat_notlog[g,CTRL4_interface]),
                  mean(ctrl1$mat_notlog[g,CTRL1_medulla]), mean(ctrl2$mat_notlog[g,CTRL2_medulla]), 
                  mean(ctrl3$mat_notlog[g,CTRL3_medulla]), mean(ctrl4$mat_notlog[g,CTRL4_medulla]) )

df<-rbind(df1,df2,df3, df4)
head(df)
df$group<-factor(df$group, levels=c("Global","Cortex","Interface","Medulla"))

ctrl_gline <- ggplot(df, aes(x = t, y = gexp, col=group)) + geom_point() + 
  scale_color_manual(values=c("Global"="black","Cortex"="red","Interface"="blue","Medulla"="green")) +
  ylim(overall_min0, overall_max0) + theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank() ) +
  xlab("Kidney_ID") + ylab("Avg. gene exp. (counts)")

### AKI

overall_min <- min(AKI_sham$mat_notlog[g,], AKI_4h$mat_notlog[g,], 
                   AKI_12h$mat_notlog[g,], AKI_2d$mat_notlog[g,])

overall_max <- quantile(c(AKI_sham$mat_notlog[g,], AKI_4h$mat_notlog[g,], 
                          AKI_12h$mat_notlog[g,], AKI_2d$mat_notlog[g,]), 0.95)

color_scale <- scale_color_gradientn(colors=c("darkgray","yellow","red"),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1), 
                                     limits = c(overall_min, overall_max), oob = scales::squish)

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)

my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

df <- data.frame(AKI_sham$pos, gexp=AKI_sham$mat_notlog[g,])
g1 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('AKI sham')
df <- data.frame(AKI_4h$pos, gexp=AKI_4h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme + ggtitle('AKI 4h')
df <- data.frame(AKI_12h$pos, gexp=AKI_12h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale + my_theme  + ggtitle('AKI 12h')
df <- data.frame(AKI_2d$pos, gexp=AKI_2d$mat_notlog[g,])
g4 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale +  theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('AKI 2d')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=V5, y=V6, col=gexp)) + geom_point(size=0.6, alpha=1) + color_scale  + my_theme + ggtitle('AKI 2d')

aki_comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))


df1 <- data.frame(t = c(0, 4, 12, 48),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d])), group = 'Cortex')

df2 <- data.frame(t = c(0, 4, 12, 48),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d])), group = 'Interface')

df3 <- data.frame(t = c(0, 4, 12, 48),
                  gexp = c(mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                           mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d])), group = 'Medulla')

df4 <- data.frame(t = c(0, 4, 12, 48),
                  gexp = c(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                           mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,])), group = 'Global')

overall_min0<-min(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                  mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,]),
                  mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d]),
                  mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d]),
                  mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d]) )

overall_max0<-max(mean(AKI_sham$mat_notlog[g,]), mean(AKI_4h$mat_notlog[g,]), 
                  mean(AKI_12h$mat_notlog[g,]), mean(AKI_2d$mat_notlog[g,]),
                  mean(AKI_sham$mat_notlog[g,AKI_cortex.sham]), mean(AKI_4h$mat_notlog[g,AKI_cortex.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_cortex.12h]), mean(AKI_2d$mat_notlog[g,AKI_cortex.2d]),
                  mean(AKI_sham$mat_notlog[g,AKI_interface.sham]), mean(AKI_4h$mat_notlog[g,AKI_interface.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_interface.12h]), mean(AKI_2d$mat_notlog[g,AKI_interface.2d]),
                  mean(AKI_sham$mat_notlog[g,AKI_medulla.sham]), mean(AKI_4h$mat_notlog[g,AKI_medulla.4h]), 
                  mean(AKI_12h$mat_notlog[g,AKI_medulla.12h]), mean(AKI_2d$mat_notlog[g,AKI_medulla.2d]) )

df<-rbind(df1,df2,df3, df4)
head(df)
df$group<-factor(df$group, levels=c("Global","Cortex","Interface","Medulla"))

aki_gline <- ggplot(df, aes(x = t, y = gexp, col=group)) + geom_line() + geom_point() + 
  scale_color_manual(values=c("Global"="black","Cortex"="red","Interface"="blue","Medulla"="green")) +
  ylim(overall_min0, overall_max0) + theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank() ) +
  xlab("time (hours)") + ylab("Avg. gene exp. (counts)")


## CIS

overall_min <- min(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                   CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,])

overall_max <- quantile(c(CIS_0h$mat_notlog[g,], CIS_12h$mat_notlog[g,], 
                          CIS_24h$mat_notlog[g,], CIS_48h$mat_notlog[g,]), 0.95)

color_scale <- scale_color_gradientn(colors=c("darkgray","yellow","red"),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1), 
                                     limits = c(overall_min, overall_max), oob = scales::squish)

color_scale <- scale_color_gradientn(colors=viridis(5),
                                     
                                     #colors=c("grey","magenta","yellow","red"),
                                     values = c(0, 0.8, 1),
                                     limits = c(overall_min, overall_max), oob = scales::squish)


my_theme <- theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(),
        legend.position="none"
  )

df <- data.frame(CIS_0h$pos, gexp=CIS_0h$mat_notlog[g,])
g1 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme  + ggtitle('CIS 0h') 
df <- data.frame(CIS_12h$pos, gexp=CIS_12h$mat_notlog[g,])
g2 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1)  + color_scale  + my_theme + ggtitle('CIS 12h')
df <- data.frame(CIS_24h$pos, gexp=CIS_24h$mat_notlog[g,])
g3 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale + my_theme + ggtitle('CIS 24h')
df <- data.frame(CIS_48h$pos, gexp=CIS_48h$mat_notlog[g,])
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + theme(legend.text = element_text(size=10), legend.ticks = element_line(color = "black")) + ggtitle('CIS 48h')

# Extract the legend from one plot
legend <- get_legend(g4)
g4 <- ggplot(df, aes(x=X, y=Y, col=gexp)) + geom_point(size=0.5, alpha=1) + color_scale  + my_theme + ggtitle('CIS 48h')

cis_comb_plot<-plot_grid(g1, g2, g3, g4, legend, ncol = 5, rel_widths = c(1, 1, 1, 1, 0.25))

df1 <- data.frame(t = c(0, 12, 24, 48),
                  gexp = c(mean(CIS_0h$mat_notlog[g,CIS_cortex.0h]), mean(CIS_12h$mat_notlog[g,CIS_cortex.12h]), 
                           mean(CIS_24h$mat_notlog[g,CIS_cortex.24h]), mean(CIS_48h$mat_notlog[g,CIS_cortex.48h])), group = 'Cortex')
df2 <- data.frame(t = c(0, 12, 24, 48),
                  gexp = c(mean(CIS_0h$mat_notlog[g,CIS_interface.0h]), mean(CIS_12h$mat_notlog[g,CIS_interface.12h]), 
                           mean(CIS_24h$mat_notlog[g,CIS_interface.24h]), mean(CIS_48h$mat_notlog[g,CIS_interface.48h])), group = 'Interface')

df3 <- data.frame(t = c(0, 12, 24, 48),
                  gexp = c(mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                           mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h])), group = 'Medulla')

df4 <- data.frame(t = c(0, 12, 24, 48),
                  gexp = c(mean(CIS_0h$mat_notlog[g,]), mean(CIS_12h$mat_notlog[g,]), 
                           mean(CIS_24h$mat_notlog[g,]), mean(CIS_48h$mat_notlog[g,])), group = 'Global')

overall_min0<-min(mean(CIS_0h$mat_notlog[g,]), mean(CIS_12h$mat_notlog[g,]), 
                  mean(CIS_24h$mat_notlog[g,]), mean(CIS_48h$mat_notlog[g,]),
                  mean(CIS_0h$mat_notlog[g,CIS_cortex.0h]), mean(CIS_12h$mat_notlog[g,CIS_cortex.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_cortex.24h]), mean(CIS_48h$mat_notlog[g,CIS_cortex.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_interface.0h]), mean(CIS_12h$mat_notlog[g,CIS_interface.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_interface.24h]), mean(CIS_48h$mat_notlog[g,CIS_interface.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h]) )

overall_max0<-max(mean(CIS_0h$mat_notlog[g,]), mean(CIS_12h$mat_notlog[g,]), 
                  mean(CIS_24h$mat_notlog[g,]), mean(CIS_48h$mat_notlog[g,]),
                  mean(CIS_0h$mat_notlog[g,CIS_cortex.0h]), mean(CIS_12h$mat_notlog[g,CIS_cortex.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_cortex.24h]), mean(CIS_48h$mat_notlog[g,CIS_cortex.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_interface.0h]), mean(CIS_12h$mat_notlog[g,CIS_interface.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_interface.24h]), mean(CIS_48h$mat_notlog[g,CIS_interface.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h]),
                  mean(CIS_0h$mat_notlog[g,CIS_medulla.0h]), mean(CIS_12h$mat_notlog[g,CIS_medulla.12h]), 
                  mean(CIS_24h$mat_notlog[g,CIS_medulla.24h]), mean(CIS_48h$mat_notlog[g,CIS_medulla.48h]) )

df<-rbind(df1,df2,df3, df4)
head(df)
df$group<-factor(df$group, levels=c("Global","Cortex","Interface","Medulla"))
# gline <- ggplot(df, aes(x = t, y = gexp, col=group)) + geom_line() + geom_point() + 
#          ylim(overall_min0, overall_max0) + theme_light() +
#          theme(axis.text.x = element_text(face = "bold", color="black"),  
#                axis.text.y = element_text(face = "bold", color="black"),
#                panel.border = element_rect(color = "black", fill = NA, size = 1),
#                legend.text = element_text(face = "bold"),
#                panel.grid.minor = element_blank() ) +
#          xlab("time (hours)") + ylab("Avg. gene exp. (counts)")


cis_gline <- ggplot(df, aes(x = t, y = gexp, col=group)) + geom_line() + geom_point() + 
  scale_color_manual(values=c("Global"="black","Cortex"="red","Interface"="blue","Medulla"="green")) +
  ylim(overall_min0, overall_max0) + theme_light() +
  theme(axis.text.x = element_text(color="black", size=11),  
        axis.text.y = element_text(color="black",size=11),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(color="black", size=12),
        panel.grid.minor = element_blank() ) +
  xlab("time (hours)") + ylab("Avg. gene exp. (counts)")


library(gridExtra)
grid.arrange(ctrl_comb_plot,ctrl_gline,
             aki_comb_plot,aki_gline,
             cis_comb_plot,cis_gline,
             nrow=3,ncol=2, widths=c(16, 5), top=g)