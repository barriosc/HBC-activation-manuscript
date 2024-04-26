##call libraries necessary 
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(grid)
library(dplyr)
source("https://goo.gl/4mthoF")
library(tidyverse)
library(DescTools)


#### Read in PMA dataset ###
masterdata <- read.csv("PMATimeCourseMaster_Data.csv", header=TRUE)

### Normalize P63 and Rela to the area of the nucleus and cytoplasm

FirstData <- transform(masterdata, NormNucP63 =  NucleusIntensity_IntegratedIntensity_P63  / NucleusAreaShape_Area)
SecondData <- transform(FirstData, NormNucRela =  NucleusIntensity_IntegratedIntensity_RELA  / NucleusAreaShape_Area)
ThirdData <- transform(SecondData, NormCytoP63=  Cytoplasm.Intensity_IntegratedIntensity_P63  / Cytoplasm.AreaShape_Area)
FinalData <- transform(ThirdData, NormCytoRelA =  Cytoplasm.Intensity_IntegratedIntensity_RELA  / Cytoplasm.AreaShape_Area)

## Remove all outliers ##


### to skip this section, load in RDS file, "HBC_PMA_Activation_TimeCourse.RDS"
## Outlier Removal for NucP63 ####

Vehicle <- subset(FinalData, TimePoint=="0",
                  select=TimePoint:NormCytoRelA)

outlierKD(Vehicle, NormNucP63)

Zero.5 <- subset(FinalData, TimePoint=="0.5",
                 select=TimePoint:NormCytoRelA)

outlierKD(Zero.5, NormNucP63)

One <- subset(FinalData, TimePoint=="1",
              select=TimePoint:NormCytoRelA)

outlierKD(One, NormNucP63)


Three <- subset(FinalData, TimePoint=="3",
                select=TimePoint:NormCytoRelA)

outlierKD(Three, NormNucP63)

Six <- subset(FinalData, TimePoint=="6",
              select=TimePoint:NormCytoRelA)

outlierKD(Six, NormNucP63)

Eight <- subset(FinalData, TimePoint=="8",
                select=TimePoint:NormCytoRelA)

outlierKD(Eight, NormNucP63)

Twelve <- subset(FinalData, TimePoint=="12",
                 select=TimePoint:NormCytoRelA)

outlierKD(Twelve, NormNucP63)

FinalData <- rbind(Vehicle, Zero.5, One, Three, Six, Eight, Twelve)

##Outlier Removal for NucRelA and CytoRelA

####NucRelA###

attach(FinalData)
Vehicle_RelA <- subset(FinalData, TimePoint=="0",
                       select=TimePoint:NormCytoRelA)

outlierKD(Vehicle_RelA, NormNucRela)

Zero.5_RelA <- subset(FinalData, TimePoint=="0.5",
                      select=TimePoint:NormCytoRelA)

outlierKD(Zero.5_RelA, NormNucRela)

One_RelA <- subset(FinalData, TimePoint=="1",
                   select=TimePoint:NormCytoRelA)

outlierKD(One_RelA, NormNucRela)


Three_RelA <- subset(FinalData, TimePoint=="3",
                     select=TimePoint:NormCytoRelA)

outlierKD(Three_RelA, NormNucRela)

Six_RelA <- subset(FinalData, TimePoint=="6",
                   select=TimePoint:NormCytoRelA)

outlierKD(Six_RelA, NormNucRela)

Eight_RelA <- subset(FinalData, TimePoint=="8",
                     select=TimePoint:NormCytoRelA)

outlierKD(Eight_RelA, NormNucRela)

Twelve_RelA <- subset(FinalData, TimePoint=="12",
                      select=TimePoint:NormCytoRelA)

outlierKD(Twelve_RelA, NormNucRela)

FinalData <- rbind(Vehicle_RelA, Zero.5_RelA, One_RelA, Three_RelA, 
                   Six_RelA, Eight_RelA, Twelve_RelA)

####CytoP63####
Vehicle_CytoP63 <- subset(FinalData,TimePoint=="0",
                          select=TimePoint:NormCytoRelA)

outlierKD(Vehicle_CytoP63, NormCytoP63)

Zero.5_CytoP63 <- subset(FinalData, TimePoint=="0.5",
                         select=TimePoint:NormCytoRelA)

outlierKD(Zero.5_CytoP63, NormCytoP63)

One_CytoP63 <- subset(FinalData, TimePoint=="1",
                      select=TimePoint:NormCytoRelA)

outlierKD(One_CytoP63, NormCytoP63)


Three_CytoP63 <- subset(FinalData, TimePoint=="3",
                        select=TimePoint:NormCytoRelA)

outlierKD(Three_CytoP63, NormCytoP63)

Six_CytoP63 <- subset(FinalData, TimePoint=="6",
                      select=TimePoint:NormCytoRelA)

outlierKD(Six_CytoP63, NormCytoP63)

Eight_CytoP63 <- subset(FinalData, TimePoint=="8",
                        select=TimePoint:NormCytoRelA)

outlierKD(Eight_CytoP63, NormCytoP63)

Twelve_CytoP63 <- subset(FinalData, TimePoint=="12",
                         select=TimePoint:NormCytoRelA)

outlierKD(Twelve_CytoP63, NormCytoP63)

FinalData <- rbind(Vehicle_CytoP63, Zero.5_CytoP63, One_CytoP63, Three_CytoP63, 
                   Six_CytoP63, Eight_CytoP63, Twelve_CytoP63)

## adding NCI to final data
attach(FinalData)
FinalDatas <- transform(FinalData, NCI_P63 = ((NormNucP63) / (NormNucP63 + NormCytoP63)))
FinalDataHBC_TimeCourse <- transform(FinalDatas, NCI_RelA = ((NormNucRela) / (NormNucRela + NormCytoRelA))) 

saveRDS(FinalDataHBC_TimeCourse, "HBC_PMA_Activation_TimeCourse.rds")

## summary data for cleaned up data
## can load RDS here!

FinalDataHBC_TimeCourse <- readRDS("HBC_PMA_Activation_TimeCourse.rds")

FinalDataHBC_TimeCourse %>%
  group_by(TimePoint)  %>% count()

## Calculation of quantiles to set threshold values

options(pillar.sigfig = 7)

FinalDataHBC_TimeCourse %>% filter(NormNucP63 > 0) %>% group_by(TimePoint) %>% 
  summarize(first_25th=quantile(NormNucP63,probs=0.25),
            second_50th=quantile(NormNucP63,probs=0.5),
            third_75th=quantile(NormNucP63,probs=0.75))


### plots

## Figure 1E, degradation of TP63 over time

FINALCellProfiler_NucP63Deg_12.png <- ggplot(FinalDataHBC_TimeCourse) + 
  aes(x = as.factor(TimePoint), y = NormNucP63, fill=as.factor(TimePoint)) + 
  geom_hline(yintercept = c(0.3397729, 0.2740572), linetype = 'dashed', color = 'black') +
  geom_violin(scale = 'width', color=NA) + theme_classic(base_line_size =1, base_rect_size = 1) +
  geom_boxplot(outlier.shape = NA, width=0.2, fill='white')+
  scale_fill_viridis_d() + xlab("Hours post-treatment (HPT)") + ylab(expression(paste(frac(Nuclear~P63~Integrated~Density, Nuclear~Area))))+
  theme(legend.position = 'none') +
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15))

final_PMA_timecourse <- FINALCellProfiler_NucP63Deg_12.png + theme(axis.text.x = element_text(color="black", 
                                                                                             size=20),
                                                                  axis.text.y = element_text(color="black", 
                                                                                             size=23), 
                                                                  axis.title.x=element_text(color="black", 
                                                                                            size=27), 
                                                                  axis.title.y=element_text(color="black", size=27))



ggsave(filename="ActivationPaper_Figure1E.png", plot=final_PMA_timecourse, device="png",
       height=6, width=10, units="in", dpi=1000)



### Dunnett Test to determine significance 

DunnettTest(x=FinalDataHBC_TimeCourse$NormNucP63, g=FinalDataHBC_TimeCourse$TimePoint)



## Granularity Data (Figure 1F)

### granularity data ###

nucleargranularitydata <- read.csv("ActivationPaper_GranularityNuclearMasks.csv", header=TRUE)  

Vehicle_ <- subset(nucleargranularitydata, TimePoint=="0",
                   select=TimePoint:Number_Object_Number)

outlierKD(Vehicle_, Granularity_1_P63)

Zero.5_ <- subset(nucleargranularitydata, TimePoint=="0.5",
                  select=TimePoint:Number_Object_Number)

outlierKD(Zero.5_, Granularity_1_P63)

One_ <- subset(nucleargranularitydata, TimePoint=="1",
               select=TimePoint:Number_Object_Number)

outlierKD(One_, Granularity_1_P63)


Three_ <- subset(nucleargranularitydata, TimePoint=="3",
                 select=TimePoint:Number_Object_Number)

outlierKD(Three_, Granularity_1_P63)

Six_ <- subset(nucleargranularitydata, TimePoint=="6",
               select=TimePoint:Number_Object_Number)

outlierKD(Six_, Granularity_1_P63)

Eight_ <- subset(nucleargranularitydata, TimePoint=="8",
                 select=TimePoint:Number_Object_Number)


outlierKD(Eight_, Granularity_1_P63)

Twelve_ <- subset(nucleargranularitydata, TimePoint=="12",
                  select=TimePoint:Number_Object_Number)

outlierKD(Twelve_, Granularity_1_P63)

nucleargranularitydata_final <- rbind(Vehicle_, Zero.5_, One_, Three_, Six_, Eight_, Twelve_)

granularity_df <- as.data.frame(nucleargranularitydata_final)

### Calculate mean and standard deviation for barplot

df_mean_std_granularity <- nucleargranularitydata_final %>%
  group_by(TimePoint) %>% filter(Granularity_1_P63 > 0) %>%
  summarise_at(vars(Granularity_1_P63), list(mean=mean, sd=sd)) %>% 
  as.data.frame()

### data visualization

p <- ggplot(df_mean_std_granularity, aes(factor(TimePoint), y = mean, fill=as.factor(TimePoint))) + 
  geom_bar(stat="identity")

p2 <- p + theme_classic() + scale_fill_viridis(discrete=TRUE, option="D")

p3 <- p2 + xlab("Hours post-treatment (HPT)") + ylab("Granularity Index") + theme_classic(base_line_size =1, base_rect_size = 1)

final <- p3 +   theme(axis.title = element_text(size=15), axis.text = element_text(size=15)) +  theme(legend.position = 'none') + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3)

final_granularity_barplot <- final + theme(axis.text.x = element_text(color="black", 
                                                       size=20),
                            axis.text.y = element_text(color="black", 
                                                       size=23), 
                            axis.title.x=element_text(color="black", 
                                                      size=20), 
                            axis.title.y=element_text(color="black", size=23))


## counts per group 

df_mean_std <- nucleargranularitydata_final %>%
  group_by(TimePoint)  %>% count()


ggsave(filename="ActivationPaper_Granularity_Figure1F.png", plot=final_granularity_barplot, device="png",
       height=4, width=6, units="in", dpi=1000)


### Figure 1G ####

## Determine the proportion of cells that belong to fluorescence-classified cells

# set the same fluorescent values as established above
FinalDataHBC_TimeCourse <- FinalDataHBC_TimeCourse %>% mutate(p63_category = 'resting')
FinalDataHBC_TimeCourse$p63_category[FinalDataHBC_TimeCourse$NormNucP63 < 0.3397729] <- 'intermediate'
FinalDataHBC_TimeCourse$p63_category[FinalDataHBC_TimeCourse$NormNucP63 < 0.2544129] <- 'activated'


PMA_treated_HBCs_categories <- FinalDataHBC_TimeCourse %>%
  group_by(TimePoint, p63_category) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

write.csv(PMA_treated_HBCs_categories, "PMA_treated_HBCs_categories.csv")

### manually added 0 to for 12HPT, intermediate category

plot <- ggplot(PMA_treated_HBCs_categories, aes(x = as.factor(TimePoint), y = percentage, fill = p63_category)) +
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


ggsave(filename="ActivationPaper_1G.png", plot=plot, device="png",
       height=5, width=5, units="in", dpi=1000)


### Figure 1H, PMA recovery assay ###

nucleardata <- read.csv("MyExpt_ActivationPaper_PMARecoveryNuclearMasks (1).csv", header=TRUE)

### Normalize P63 to the area of the nucleus and cytoplasm

NormNucP63 <- transform(nucleardata, NormNucP63 =  Intensity_IntegratedIntensity_P63  / AreaShape_Area)

## Remove all outliers ##

## Outlier Removal for NucP63 ####

Vehicle <- subset(NormNucP63, TimePoint=="0",
                  select=TimePoint:NormNucP63)

outlierKD(Vehicle, NormNucP63)

Six <- subset(NormNucP63, TimePoint=="6",
              select=TimePoint:NormNucP63)

outlierKD(Six, NormNucP63)

Twelve <- subset(NormNucP63, TimePoint=="12",
                 select=TimePoint:NormNucP63)

outlierKD(Twelve, NormNucP63)
NormNucP63 <- rbind(Vehicle, Six, Twelve)

DunnettTest(x=NormNucP63$NormNucP63, g=NormNucP63$TimePoint)

#define quantiles of interest
options(pillar.sigfig = 7)
q = c(.25, .5, .75)

NormNucP63 %>% filter(NormNucP63 > 0) %>% group_by(TimePoint) %>% 
  summarize(quant25 = quantile(NormNucP63, probs = q[1]), 
            quant50 = quantile(NormNucP63, probs = q[2]),
            quant75 = quantile(NormNucP63, probs = q[3]))

NormNucP63 <- NormNucP63 %>% mutate(p63_category = 'resting')
NormNucP63$p63_category[NormNucP63$NormNucP63 < 0.3338766 & NormNucP63$NormNucP63 > 0.2716799 ] <- 'intermediate'
NormNucP63$p63_category[NormNucP63$NormNucP63 < 0.2716799] <- 'activated'


final_recovery <- ggplot(NormNucP63) + 
  aes(x = as.factor(TimePoint), y = NormNucP63, fill=as.factor(TimePoint)) + 
  geom_hline(yintercept = c(0.3338766, 0.2716799), linetype = 'dashed', color = 'black') +
  geom_violin(scale = 'width', color=NA) + theme_classic(base_line_size =1, base_rect_size = 1) +
  geom_boxplot(outlier.shape = NA, width=0.2, fill='white')+
  scale_fill_viridis_d() + xlab("Treatment time (hours) prior to recovery (TPR)") + 
  ylab(expression(paste(frac(Nuclear~P63~Integrated~Density, Nuclear~Area)))) +
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(color="black", 
                                   size=20),
        axis.text.y = element_text(color="black", 
                                   size=23), 
        axis.title.x=element_text(color="black", 
                                  size=20), 
        axis.title.y=element_text(color="black", size=23))



ggsave(filename="ActivationPaper_Figure1I.png", plot=final_recovery, device="png",
       height=6, width=8, units="in", dpi=1000)


## Figure IJ ## 


PMA_recovery_counts_category <- NormNucP63 %>%
  group_by(TimePoint, p63_category) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)



plot <- ggplot(PMA_recovery_counts_category, aes(x = as.factor(TimePoint), y = percentage, fill = p63_category)) +
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


ggsave(filename="ActivationPaper_Figure1J.png", plot=plot, device="png",
       height=5, width=5, units="in", dpi=1000)



