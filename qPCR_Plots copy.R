#Clean environment
rm(list=ls(all.names=TRUE))
gc()

library(ggplot2)
library(ggsci)

workdir <- setwd("/Users/srujansingh/Library/CloudStorage/OneDrive-SharedLibraries-JohnsHopkins/Jean Fan - 2024_spatial_AKI")

g <- 'Atp5b'
data <- read.csv(paste0(workdir, "/Ryo_Data/", g, ".csv"))
#data$Time<-gsub("h","",data$Time)
data$Region <- factor(data$Region, levels=c("Cortex", "Medulla"))
data$Time <- factor(data$Time, levels=c("0h", "12h","24h","48h"))


## Ryo's plot
graph1 <- ggplot(data, aes(x=Time, y=RQ, color=Region,group=Region, fill=Region))+
  stat_summary(fun.data="mean_se",geom="errorbar", position=position_dodge(width=0.9),width=0.4, show.legend=F)+
  stat_summary(fun="mean", geom="bar", position="dodge",show.legend=F,alpha=0.2)+
  geom_point(position=position_jitterdodge(jitter.width=0.1,jitter.height=0,dodge.width=0.9),alpha=1)+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,2))+
  scale_fill_aaas(name="",labels=c("Cortex","Medulla"))+
  scale_color_aaas(name="")+
  labs(y="Transcription level", title ="Atp5b")+
  theme_classic()+
  theme(plot.title=element_text(size=9,face="bold.italic", hjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=8)
  )+
  geom_text(aes(x=2,y=1.9,label="*"), color="black")+
  geom_segment(aes(x=1.7,xend=2.3,y=1.85,yend=1.85),color="black")

graph1

## sami did log (hours+1)
#data$Time2 <- data$Time
#levels(data$Time2) = c("0", "12", "24", "48")
#data$Time2 <- as.numeric(data$Time2)
#data$Time2 <- log10(data$Time2 + 1)
## doesn't make much of a difference in visualization

## Jean facet
ggplot(data, aes(x=Time, y=RQ, color=Region, group=Region, fill=Region)) + 
  #geom_point(position=position_jitterdodge(jitter.width=0.1,jitter.height=0,dodge.width=0.9),alpha=1)+
  geom_point() + 
  #geom_violin() + 
  #scale_y_continuous(expand=c(0,0))+
  #coord_cartesian(ylim=c(0,2))+
  geom_smooth(method='lm') + 
  scale_fill_aaas(name="",labels=c("Cortex","Medulla"))+
  scale_color_aaas(name="")+
  labs(y="Transcription level", title =g)+
  theme_classic()+
  theme(plot.title=element_text(size=9,face="bold.italic", hjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=8)
  ) + facet_wrap(nrow=2, vars(Region), scales="free_y")

#Plotted on the same plot
ggplot(data, aes(x = Time, y = RQ, color = Region, group = Region, fill = Region)) +  
  geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
  geom_smooth(method = 'lm', se = TRUE) +  
  scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
  scale_color_aaas(name = "") +  
  labs(y = "Transcription level", title = g) +  
  theme_classic() +  
  theme(plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5),  
        axis.title.x = element_blank(),  
        axis.title.y = element_text(size = 8),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),  
        axis.ticks.x = element_blank(),  
        legend.text = element_text(size = 8))


#### Cox6a1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g <- 'Cox6a'
data <- read.csv(paste0(workdir, "/Ryo_Data/", g, ".csv"))
#data$Time<-gsub("h","",data$Time)
data$Region <- factor(data$Region, levels=c("Cortex", "Medulla"))
data$Time <- factor(data$Time, levels=c("0h", "12h","24h","48h"))


#Plotted on the same plot
ggplot(data, aes(x = Time, y = RQ, color = Region, group = Region, fill = Region)) +  
  geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
  geom_smooth(method = 'lm', se = TRUE) +  
  scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
  scale_color_aaas(name = "") +  
  labs(y = "Transcription level", title = g) +  
  theme_classic() +  
  theme(plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5),  
        axis.title.x = element_blank(),  
        axis.title.y = element_text(size = 8),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),  
        axis.ticks.x = element_blank(),  
        legend.text = element_text(size = 8))


#### Ppargc1a ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g <- 'Pgc1a'
data <- read.csv(paste0(workdir, "/Ryo_Data/", g, ".csv"))
#data$Time<-gsub("h","",data$Time)
data$Region <- factor(data$Region, levels=c("Cortex", "Medulla"))
data$Time <- factor(data$Time, levels=c("0h", "12h","24h","48h"))


#Plotted on the same plot
ggplot(data, aes(x = Time, y = RQ, color = Region, group = Region, fill = Region)) +  
  geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
  geom_smooth(method = 'lm', se = TRUE) +  
  scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
  scale_color_aaas(name = "") +  
  labs(y = "Transcription level", title = g) +  
  theme_classic() +  
  theme(plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5),  
        axis.title.x = element_blank(),  
        axis.title.y = element_text(size = 8),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),  
        axis.ticks.x = element_blank(),  
        legend.text = element_text(size = 8))


#### Pink1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g <- 'Pink1'
data <- read.csv(paste0(workdir, "/Ryo_Data/", g, ".csv"))
#data$Time<-gsub("h","",data$Time)
data$Region <- factor(data$Region, levels=c("Cortex", "Medulla"))
data$Time <- factor(data$Time, levels=c("0h", "12h","24h","48h"))


#Plotted on the same plot
ggplot(data, aes(x = Time, y = RQ, color = Region, group = Region, fill = Region)) +  
  geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
  geom_smooth(method = 'lm', se = TRUE) +  
  scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
  scale_color_aaas(name = "") +  
  labs(y = "Transcription level", title = g) +  
  theme_classic() +  
  theme(plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5),  
        axis.title.x = element_blank(),  
        axis.title.y = element_text(size = 8),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),  
        axis.ticks.x = element_blank(),  
        legend.text = element_text(size = 8))


#### Irf3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g <- 'Irf3'
data <- read.csv(paste0(workdir, "/Ryo_Data/", g, ".csv"))
#data$Time<-gsub("h","",data$Time)
data$Region <- factor(data$Region, levels=c("Cortex", "Medulla"))
data$Time <- factor(data$Time, levels=c("0h", "12h","24h","48h"))


#Plotted on the same plot
ggplot(data, aes(x = Time, y = RQ, color = Region, group = Region, fill = Region)) +  
  geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
  geom_smooth(method = 'lm', se = TRUE) +  
  scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
  scale_color_aaas(name = "") +  
  labs(y = "Transcription level", title = g) +  
  theme_classic() +  
  theme(plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5),  
        axis.title.x = element_blank(),  
        axis.title.y = element_text(size = 8),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),  
        axis.ticks.x = element_blank(),  
        legend.text = element_text(size = 8))


#### Casp3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g <- 'Casp3'
data <- read.csv(paste0(workdir, "/Ryo_Data/", g, ".csv"))
#data$Time<-gsub("h","",data$Time)
data$Region <- factor(data$Region, levels=c("Cortex", "Medulla"))
data$Time <- factor(data$Time, levels=c("0h", "12h","24h","48h"))


#Plotted on the same plot
ggplot(data, aes(x = Time, y = RQ, color = Region, group = Region, fill = Region)) +  
  geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
  geom_smooth(method = 'lm', se = TRUE) +  
  scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
  scale_color_aaas(name = "") +  
  labs(y = "Transcription level", title = g) +  
  theme_classic() +  
  theme(plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5),  
        axis.title.x = element_blank(),  
        axis.title.y = element_text(size = 8),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),  
        axis.ticks.x = element_blank(),  
        legend.text = element_text(size = 8))

#### Fth1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g <- 'Fth1'
data <- read.csv(paste0(workdir, "/Ryo_Data/", g, ".csv"))
#data$Time<-gsub("h","",data$Time)
data$Region <- factor(data$Region, levels=c("Cortex", "Medulla"))
data$Time <- factor(data$Time, levels=c("0h", "12h","24h","48h"))


#Plotted on the same plot
ggplot(data, aes(x = Time, y = RQ, color = Region, group = Region, fill = Region)) +  
  geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
  geom_smooth(method = 'lm', se = TRUE) +  
  scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
  scale_color_aaas(name = "") +  
  labs(y = "Transcription level", title = g) +  
  theme_classic() +  
  theme(plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5),  
        axis.title.x = element_blank(),  
        axis.title.y = element_text(size = 8),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),  
        axis.ticks.x = element_blank(),  
        legend.text = element_text(size = 8))


#############################

# library(ggplot2)
# library(ggpubr)
# 
# # Perform pairwise comparisons between Cortex and Medulla at each Time point
# stat_test <- data %>%
#   group_by(Time) %>%
#   rstatix::t_test(RQ ~ Region) %>%
#   rstatix::adjust_pvalue(method = "bonferroni") %>%
#   rstatix::add_significance("p.adj") %>% # Ensure it's a data frame before using mutate()
#   mutate(y.position = max((data.frame(data)$RQ), na.rm = TRUE) * 1.1)  # Add y position  # Place p-value above max y-value
# 
# 
# ggplot(data, aes(x = factor(Time, levels = c("0", "12", "24", "48")),  
#                  y = RQ, color = factor(Region, levels = c("Cortex", "Medulla")),  
#                  group = factor(Region, levels = c("Cortex", "Medulla")),  
#                  fill = factor(Region, levels = c("Cortex", "Medulla")))) +  
#   geom_point(position = position_dodge(width = 0.2), alpha = 1) +  
#   geom_smooth(method = 'lm', se = TRUE) +  
#   stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.2), 
#                shape = 21, size = 3, stroke = 1.2) +  # Show mean points
#   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, 
#                position = position_dodge(width = 0.2)) +  # Show error bars
#   scale_fill_aaas(name = "", labels = c("Cortex", "Medulla")) +  
#   scale_color_aaas(name = "") +  
#   labs(y = "Transcription level", title = g) +  
#   theme_classic() +  
#   theme(plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5),  
#         axis.title.x = element_blank(),  
#         axis.title.y = element_text(size = 8),  
#         axis.text.x = element_text(size = 7),  
#         axis.text.y = element_text(size = 7),  
#         axis.ticks.x = element_blank(),  
#         legend.text = element_text(size = 8)) +
#   stat_pvalue_manual(data.frame(stat_test), label = "p.adj.signif",  
#                      tip.length = 0.02,  
#                      bracket.size = 0.5, xmin=stat_test$group1, xmax=stat_test$group2)
# 
