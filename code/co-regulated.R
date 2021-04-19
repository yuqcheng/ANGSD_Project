suppressPackageStartupMessages({
  library(DESeq2)
  library(magrittr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(ggpubr)
  library(ggthemes)
  library(eulerr)
})

DESeq.ds <- readRDS('~/R Scripts/angsd-project/featurecount_ds.rds')
dds <- DESeq(DESeq.ds)

#NRG1_regulated
res_NRG1_DMSO <- results(dds, contrast = c("condition","DMSO", "NRG1"))
res_NRG1_DMSO$Symbol <- rownames(res_NRG1_DMSO)
#AR¡ª¡ªregulated
res_Enz_DMSO <- results(dds, contrast = c("condition","DMSO", "Enz"))
res_Enz_DMSO$Symbol <- rownames(res_Enz_DMSO)

res_Enz_DMSO <- subset(res_Enz_DMSO, padj < 0.05 & abs(log2FoldChange) > 0)
res_NRG1_DMSO <- subset(res_NRG1_DMSO, padj < 0.05 & abs(log2FoldChange) > 0)

ar <- dim(res_Enz_DMSO)
nrg <- dim(res_NRG1_DMSO)

corg <- length(which(res_Enz_DMSO$Symbol %in% res_NRG1_DMSO$Symbol))
co_gene <- rownames(res_Enz_DMSO[which(res_Enz_DMSO$Symbol %in% res_NRG1_DMSO$Symbol),])

co_heatmap <- assay(dds[co_gene,])

annotation_col <- data.frame(Treatment = c("Veh","Veh","Veh","Enz","Enz","Enz","Enz*NRG1","Enz*NRG1","Enz*NRG1","NRG1","NRG1","NRG1"))
rownames(annotation_col) <- colnames(co_heatmap)
annotation_row <- data.frame(Cluster = row_cluster)
rownames(annotation_row) <- rownames(co_heatmap)
p<-pheatmap(co_heatmap, scale = 'row', color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)), 
         legend_labels='Z-score', cluster_cols=F, 
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         angle_col = 45,show_rownames = F, cutree_rows = 4)

row_cluster <- cutree(p$tree_row,k=4)
newOrder <- co_heatmap[p$tree_row$order,]
a <- row_cluster[match(rownames(newOrder),names(row_cluster))]
length(which(a == 1))
length(which(a == 2))
length(which(a == 3))
length(which(a == 4))
colnames(newOrder)[ncol(newOrder)]="Cluster"

#co-regulated
vd <- euler(c(AR_regulated = ar[1]-corg, NRG1_regulated = nrg[1]-corg, 'AR_regulated&NRG1_regulated' = corg))
plot(vd,
     fills = list(fill = c("#ffff99", "#beaed4"), alpha = 0.7), 
     legend = list(side = "right",fontsize=15),
     edges = FALSE,
     quantities = list(fontsize = 18))

NRG1_up <- subset(res_NRG1_DMSO, log2FoldChange > 1)
NRG1_down <- subset(res_NRG1_DMSO, log2FoldChange < 1)
AR_up <- subset(res_Enz_DMSO, log2FoldChange < 1)
AR_down <- subset(res_Enz_DMSO, log2FoldChange > 1)

