### This script details the analysis done for Figure 2A, 2B, 2D, 2F, 2G, 2H, and 2I

### Single cell analyses ###

library(Seurat)

## Figure 2A

# load Gadye HBC seurat object, plot clusters
seurat_object <- readRDS('SeuratObject_Gadye_HBCs.rds')

### Figure 2A
DimPlot(seurat_object) & 
  theme_linedraw(base_line_size = 0, base_rect_size = 1) & 
  theme(axis.text = element_blank())  &
  ylab('UMAP2') & xlab('UMAP1') &
  theme(legend.position='none')

### Figure 2B

HBC_markers <- FeaturePlot(seurat_object, features = c('Krt5','Krt14','Trp63', 'Mki67', 'Hopx', "Hbegf","Krt16","Dmkn","Sbsn","Sprr1a"), 
             ncol = 2, cols = c('grey90','navyblue')) &
  theme_linedraw(base_line_size = 0, base_rect_size = 1) & 
  theme(axis.text = element_blank())  &
  ylab('UMAP2') & xlab('UMAP1')

### Figure 2D

# plot signature module expressions
FeaturePlot(seurat_object, reduction = 'umap', 
            min.cutoff = c('q5','q5'),
            max.cutoff = c('q90','q90'),
            features = c('woundresponse1','pma1'), 
            ncol=2) & 
  scale_color_gradientn(colors = c('cyan3','grey90','red2')) &
  theme_linedraw(base_line_size = 0, base_rect_size = 1) & 
  theme(legend.position = 'none') & theme(axis.text = element_blank()) &
  ylab('UMAP2') & xlab('UMAP1')

### Figures 2F and 2G

library("clusterProfiler")
library("wordcloud")
library("org.Rn.eg.db")
library("ggplot2")
organism = "org.Rn.eg.db"
library("enrichplot")
library("org.Mm.eg.db")
organism_mouse="org.Mm.eg.db"
library("cowplot")

DEGs_6HPT <- read.csv("activation_6hpt_subsetted.csv", header=TRUE)
DEGs_12HPT <- read.csv("activation_12hpt_subsetted.csv", header=TRUE)
DEGs_activation <- read.csv("Gadye_ActivationHBCs.csv", header=TRUE)
DEGs_restingHBC <-read.csv("restingHBC_gadye.csv", header=TRUE)
DEGs_HBCstar2 <- read.csv("Gadye_HBCstar2.csv", header=T)
DEGs_HBCstar1 <- read.csv("Gadye_HBCstar1.csv", header=T)

##upregulated genes and downregulated

DEGs_6HPT_up <- subset(DEGs_6HPT, log2FoldChange > 0)
DEGs_6HPT_dowm <- subset(DEGs_6HPT, log2FoldChange < 0)

DEGs_12HPT_up <- subset(DEGs_12HPT, log2FoldChange > 0)
DEGs_12HPT_dowm <- subset(DEGs_12HPT, log2FoldChange < 0)

DEGs_activation_up <- subset(DEGs_activation, logFC > 0)
DEGs_activation_dowm <- subset(DEGs_activation, logFC < 0)

DEGs_resting_up <- subset(DEGs_restingHBC, logFC > 0)
DEGs_resting_dowm <- subset(DEGs_restingHBC, logFC < 0)

DEGs_HBCstar1g_up <- subset(DEGs_HBCstar2, logFC > 0)
DEGs_HBCstar1g_dowm <- subset(DEGs_HBCstar2, logFC < 0)

DEGs_HBCstar2_up <- subset(DEGs_HBCstar1, logFC > 0)
DEGs_HBCstar2_dowm <- subset(DEGs_HBCstar1, logFC < 0)

## upcluster

DEGs_6HPT_up <- DEGs_6HPT_up$X
DEGs_12HPT_up <- DEGs_12HPT_up$X
DEGs_activation_up <- DEGs_activation_up$Feature
DEGs_resting_up <- DEGs_resting_up$Feature
DEGs_HBCstar1g_up <- DEGs_HBCstar1g_up$Feature
DEGs_HBCstar2_up <- DEGs_HBCstar2_up$Feature

DEGs_6HPT_up <- bitr(DEGs_6HPT_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_12HPT_up <- bitr(DEGs_12HPT_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_activation_up <- bitr(DEGs_activation_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_resting_up <- bitr(DEGs_resting_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_HBCstar1g_up <- bitr(DEGs_HBCstar1g_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_HBCstar2_up <- bitr(DEGs_HBCstar2_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

up_data <- list(DEGs_6HPT_upm=as.vector(unlist(DEGs_6HPT_up)), 
                DEGs_12HPT_up=as.vector(unlist(DEGs_12HPT_up)),
                DEGs_resting_up=as.vector(unlist(DEGs_resting_up)),
                DEGs_HBCstar1g_up=as.vector(unlist(DEGs_HBCstar1g_up)),
                DEGs_HBCstar2_up=as.vector(unlist(DEGs_HBCstar2_up)),
                DEGs_activation_up=as.vector(unlist(DEGs_activation_up)))


ck_up <- compareCluster(geneClusters = up_data, fun = "enrichGO",
                        OrgDb = organism_mouse, ont= "BP", pvalueCutoff  = .01)


simplified_ck_up <- simplify(
  ck_up,
  cutoff = 0.70,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)

HBC_deg_comparison_up <- dotplot(simplified_ck_up, showCategory=12, font.size = 8, label_format = 120)


ggsave(filename="ActivationPaper_Figure2F.png", plot=HBC_deg_comparison_up, device="png",
       height=6, width=10, units="in", dpi=1000)

## down

DEGs_6HPT_dowm <- DEGs_6HPT_dowm$X
DEGs_12HPT_dowm <- DEGs_12HPT_dowm$X
DEGs_activation_dowm <- DEGs_activation_dowm$Feature
DEGs_resting_dowm <- DEGs_resting_dowm$Feature
DEGs_HBCstar1g_dowm <- DEGs_HBCstar1g_dowm$Feature
DEGs_HBCstar2_dowm <- DEGs_HBCstar2_dowm$Feature

DEGs_6HPT_dowm <- bitr(DEGs_6HPT_dowm, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_12HPT_dowm <- bitr(DEGs_12HPT_dowm, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_activation_dowm <- bitr(DEGs_activation_dowm, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_resting_dowm <- bitr(DEGs_resting_dowm, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_HBCstar1g_dowm <- bitr(DEGs_HBCstar1g_dowm, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
DEGs_HBCstar2_dowm <- bitr(DEGs_HBCstar2_dowm, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

dowm_data <- list(DEGs_6HPT_dowm=as.vector(unlist(DEGs_6HPT_dowm)), 
                  DEGs_12HPT_dowm=as.vector(unlist(DEGs_12HPT_dowm)),
                  DEGs_resting_dowm=as.vector(unlist(DEGs_resting_dowm)),
                  DEGs_HBCstar1g_dowm=as.vector(unlist(DEGs_HBCstar1g_dowm)),
                  DEGs_HBCstar2_dowm=as.vector(unlist(DEGs_HBCstar2_dowm)))
#DEGs_activation_dowm=as.vector(unlist(DEGs_activation_dowm)))

ck_down <- compareCluster(geneClusters = dowm_data, fun = "enrichGO",
                          OrgDb = organism_mouse, ont= "BP", pvalueCutoff  = .01)


simplified_ck <- simplify(
  ck_down,
  cutoff = 0.70,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)

HBC_deg_comparison_down <- dotplot(ck_down, showCategory=12, font.size = 8, label_format = 120)

ggsave(filename="ActivationPaper_2G.png", plot=HBC_deg_comparison_down, device="png",
       height=6, width=10, units="in", dpi=1000)

### Figure 2H

# load in normalized expression data
normalized_counts <- readRDS('PMA_activation_DESeq2_normalized_counts.rds')

# define genes to plot
woundresponse <- c('Krt6a', 'Krt16', 'Sprr1a', 'Sprr2a3', 'Krtdap', 'Dmkn', 'Sbsn', 'Hbegf')
rat_woundresponse <- woundresponse[woundresponse %in% row.names(normalized_counts)]
heatmap_genes <- c('Tp63','Krt5','Krt14', rat_woundresponse, 'Hopx')

# define annotation dataframes
gene_ann <- data.frame(row.names = heatmap_genes, GeneSet = c( rep('Resting HBC Markers', 3),
                                                               rep('Wound Response Markers', times = length(rat_woundresponse)), 
                                                               rep('Activated HBC Marker',1)))

sample_ann <- data.frame(row.names = colnames(normalized_counts), 
                         Timepoint = paste0(rep(x =c(0,6,12), each = 3),' HPT'))
# plot heatmap
final_heatmap <- pheatmap(normalized_counts[heatmap_genes,], annotation_row = gene_ann, #annotation_col = sample_ann,
                          scale='row', cluster_cols = F, color = RColorBrewer::brewer.pal(n = 9, name = 'Greys'),
                          show_colnames = F, fontsize_row = 20, gaps_col= c(3,6), border_color = "black")

ggsave(filename="ActivationPaper_Figure2H", plot=final_heatmap, device="png",
       height=7, width=10, units="in", dpi=1000)

## Figure 2I Deconvolution





