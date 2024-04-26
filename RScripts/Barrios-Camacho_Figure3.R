## Figure 3 ##

# Figure 3A #
## read in data
total_TFs <- read.csv("rat_TFs_FINAL.csv", header=TRUE)
Vehiclevs6 <- read.csv("activation_6hpt_subsetted.csv", header=TRUE)
Vehiclevs12 <- read.csv("activation_12hpt_subsetted.csv", header=TRUE)
Latent12vs6 <- read.csv("activation_12hpt6hpt_subsetted.csv", header=TRUE)

allTFsin6HPT <- intersect(total_TFs$Symbol, Vehiclevs6$X)
allTFsin12HPT <- intersect(total_TFs$Symbol, Vehiclevs12$X)
allTFsinLate <- intersect(total_TFs$Symbol, Latent12vs6$X)

allTFs <- intersect(intersect(allTFsin6HPT,allTFsin12HPT),allTFsinLate)


## Regulon Analysis ##
library("RTN")
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(tidyr)
require(grid)

normalizedcounts <- read.csv("normalized_counts.csv", header=TRUE, row.names = 1)
matrix <- data.matrix(normalizedcounts, rownames.force = NA)

rtni <- tni.constructor(expData = matrix, 
                        regulatoryElements = allTFs)
RTNI.set <- tni.preprocess(rtni)
RTNI.set <- tni.permutation(RTNI.set, nPermutations = 100)
RTNI.set <- tni.dpi.filter(RTNI.set)
tni.regulon.summary(RTNI.set)
RTNI.set <- tni.gsea2(RTNI.set, regulatoryElements = allTFs)
metabric_regact <- tni.get(RTNI.set, what = "regulonActivity")

TF_RegulatoryProfiles <- pheatmap(t(metabric_regact$differential), 
                                  show_colnames = F, annotation_legend = TRUE, 
                                  show_rownames =F,
                                  clustering_method = "ward.D2",
                                  clustering_distance_rows = "correlation",
                                  clustering_distance_cols = "correlation",
                                  cutree_rows = 5,
                                  color= cividis(10000),
                                  cellwidth= 10,
                                  cellheight= 2,
                                  scale = "row",
                                  cluster_cols = FALSE)


ggsave(filename="TF_RegulatoryProfiles.png", plot=TF_RegulatoryProfiles, device="png",
       height=10, width=10, units="in", dpi=1000)

## GitHub has a .csv of Regulon Activity Profiles (RAP) of the matrix calculated above, GitHub also contains
## necessary cytoscape files to explore regulons and the topology analysis 

### Figure 3C ###

# missing as of 4/26/2024

## Figure 3E

## please run Figure 1 data, or read-in RDS file

## subset for 30'

FinalDatas <- readRDS("HBC_PMA_Activation_TimeCourse.rds")

subsetted <- FinalDatas %>% subset(TimePoint == "0" | TimePoint =="0.5" )

p_RelA <- ggplot(FinalDatas, aes(factor(TimePoint), y = NCI_RelA, fill=as.factor(TimePoint))) + 
  geom_violin(scale="width", trim="FALSE") + geom_boxplot(width=.2, fill="white", outlier.shape = NA, width=0.2, fill='white')

p2_RelA <- p_RelA + theme_classic() + scale_fill_viridis(discrete=TRUE, option="D")

p3_RelA <- p2_RelA + xlab("Hours post-treatment (HPT)") + ylab("RelA Nucleus-Cytoplasm Index (NCI)")


CellProfiler_RelA_NCI.png <- p3_RelA + theme(axis.text.x = element_text(color="black", 
                                                                        size=15),
                                             axis.text.y = element_text(color="black", 
                                                                        size=15), axis.title.x=element_text(color="black", 
                                                                                                            size=15), 
                                             axis.title.y=element_text(color="black", size=14))



ggsave(filename="CellProfiler_RelANuclearIntensity.png", plot=CellProfiler_RelA_Intensity.png, device="png",
       height=6, width=10, units="in", dpi=1000)


### ANOVA for nuclear translocation
summary <- aov(NCI~as.factor(TimePoint), data=FinalDatas)
TukeyHSD(summary)
set.seed(20140123)

DunnettTest(x=FinalDatas$NCI, g=FinalDatas$TimePoint)



