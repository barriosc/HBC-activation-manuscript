## The script detailed below details the analysis and visualization from Barrios-Camacho et al., Figure 1B-C

## load packages

library(RColorBrewer)
library(viridis)
library(grid)
library(plyr)
source("https://goo.gl/4mthoF")
library(tidyverse)
library(scales)

## Figure 1b AND 1C in vivo methimazole data. Please download data from RawData folder on GitHub ##

## Read in raw data

MethimazoleData <- read.csv("ActivationPaper_MethimazoleData.csv", header=TRUE)

## Calculate group quantiles to set fluorescence thresholds for HBC state ##

options(pillar.sigfig = 7)

MethimazoleData %>% filter(rescaled > 0) %>% group_by(TimePoint) %>% 
  summarize(first_25th=quantile(rescaled,probs=0.25),
            second_50th=quantile(rescaled,probs=0.5),
            third_75th=quantile(rescaled,probs=0.75))

## Data Visualization ##

p <- ggplot(MethimazoleData, aes(factor(TimePoint), y = rescaled, fill=as.factor(TimePoint))) + 
  geom_violin(scale="width", trim="FALSE") + geom_boxplot(width=.2, fill="white")

p2 <- p + theme_classic() + scale_fill_viridis(discrete=TRUE, option="D")

p3 <- p2 + xlab("Hours post-injection (HPI)") + 
  ylab(expression(paste(frac(Nuclear~P63~Integrated~Density, Nuclear~Area))))


CellProfiler_invivo_mtz_data <- p3 + theme(axis.text.x = element_text(color="black", 
                                                                      size=35),
                                           axis.text.y = element_text(color="black", 
                                                                      size=25), 
                                           axis.title.x=element_text(color="black", 
                                                                     size=35), 
                                           axis.title.y=element_text(color="black", size=35))


final_invivo <- CellProfiler_invivo_mtz_data + geom_hline(yintercept=c(0.337934 ,0.2540),linetype=2) + ylim(0,1) #+ theme_classic()

ggsave(filename="ActivationPaper_Figure1B.png", plot=final_invivo, device="png",
       height=8, width=25, units="in", dpi=1000)


## Set groups according to fluorescent values denoted above. This is the data visualization for Figure Figure 1C

MethimazoleData <- MethimazoleData %>% mutate(p63_category = 'resting')
MethimazoleData$p63_category[MethimazoleData$rescaled < 0.337934 & MethimazoleData$rescaled > 0.2540 ] <- 'transitioning'
MethimazoleData$p63_category[MethimazoleData$rescaled < 0.2540] <- 'activated'

HBC_categories_in_vivo_MTZ <- ggplot(MethimazoleData, aes(x = as.factor(TimePoint), y = rescaled, fill = p63_category)) +
  geom_col(position = "fill") +
  theme_linedraw(base_line_size = 0, base_rect_size = 1) +
  scale_fill_manual(values = rev(c('black','gold2','red2'))) +
  xlab('Timepoint') + ylab('Percent') + theme(axis.title = element_text(size=15)) + theme(legend.position = 'none')  +
  theme(axis.text.x = element_text(color="black", 
                                   size=20),
        axis.text.y = element_text(color="black", 
                                   size=23), 
        axis.title.x=element_text(color="black", 
                                  size=20), 
        axis.title.y=element_text(color="black", size=23))

ggsave(filename="ActivationPaper_Figure1C.png", plot=HBC_categories_in_vivo_MTZ, device="png",
       height=5, width=5, units="in", dpi=1000)



