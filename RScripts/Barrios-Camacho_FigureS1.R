### Supplememtal Figure S1 ###

library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(grid)
library(dplyr)
source("https://goo.gl/4mthoF")
library(tidyverse)
library(DescTools)

# read in .RDS cellprofiler data

FinalDataHBC_TimeCourse <- readRDS("HBC_PMA_Activation_TimeCourse.rds")

### For Figure S1-A

p <- ggplot(FinalDataHBC_TimeCourse, aes(factor(TimePoint), y = NCI_P63, fill=as.factor(TimePoint))) + 
  geom_violin(scale="width", trim="FALSE") + geom_boxplot(width=.2, fill="white", outlier.shape = NA, width=0.2, fill='white')

p2 <- p + theme_classic() + scale_fill_viridis(discrete=TRUE, option="D")

p3 <- p2 + xlab("Hours post-treatment (HPT)") + ylab("P63 Nuclear-Cytoplasm Index (NCI)")


NCI_P63 <- p3 + theme(axis.text.x = element_text(color="black", 
                                                                     size=20),
                                          axis.text.y = element_text(color="black", 
                                                                     size=23), axis.title.x=element_text(color="black", 
                                                                                                         size=20), 
                                          axis.title.y=element_text(color="black", size=23))

ggsave(filename="ActivationPaper_Figure1SA.png", plot=NCI_P63, device="png",
       height=6, width=10, units="in", dpi=1000)

df_mean_std <- FinalDataHBC_TimeCourse %>%
  group_by(TimePoint)  %>% count()



### For Fig S1-B
p <- ggplot(FinalDataHBC_TimeCourse, aes(factor(TimePoint), y = NormCytoP63, fill=as.factor(TimePoint))) + 
  geom_violin(scale="width", trim="FALSE") + geom_boxplot(width=.2, fill="white", outlier.shape = NA, width=0.2, fill='white')

p2 <- p + theme_classic() + scale_fill_viridis(discrete=TRUE, option="D")

p3 <- p2 + xlab("Hours post-treatment (HPT)") + ylab(expression(paste(frac(Cytoplasmic~P63~Integrated~Density, Cytoplasmic~Area))))


CellProfiler_CytoP63Deg.png <- p3 + theme(axis.text.x = element_text(color="black", 
                                                                    size=20),
                                         axis.text.y = element_text(color="black", 
                                                                    size=23), axis.title.x=element_text(color="black", 
                                                                                                        size=20), 
                                         axis.title.y=element_text(color="black", size=23))

ggsave(filename="ActivationPaper_Figure2SB.png", plot=CellProfiler_CytoP63Deg.png, device="png",
       height=6, width=10, units="in", dpi=1000)


