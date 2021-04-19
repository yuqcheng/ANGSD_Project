suppressPackageStartupMessages({
  library(DESeq2)
  library(magrittr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(ggpubr)
  library(ggthemes)
})

DESeq.ds <- readRDS('~/R Scripts/angsd-project/featurecount_ds.rds')
dds <- DESeq(DESeq.ds)
dds <- dds[rowSums(counts(dds)) > 1, ]

resultsNames(dds)

res_Enz_DMSO <- results(dds, contrast = c("condition","DMSO", "Enz"))

res_Enz_DMSO <- as.data.frame(res_Enz_DMSO)

res_Enz_DMSO$Symbol <- rownames(res_Enz_DMSO)

res_Enz_DMSO <- res_Enz_DMSO %>% dplyr::filter(abs(log2FoldChange)>0 & padj < 0.05)

res_Enz_DMSO$logP <- -log10(res_Enz_DMSO$padj)

res_Enz_DMSO$Group <- 'not-significant'

res_Enz_DMSO$Group[which(res_Enz_DMSO$log2FoldChange > 2)] <- 'Up-regulated'
res_Enz_DMSO$Group[which(res_Enz_DMSO$log2FoldChange < -2)] <- 'Down-regulated'

res_Enz_DMSO$label <- ""
res_Enz_DMSO <- res_Enz_DMSO[order(res_Enz_DMSO$padj),]
up.gene <- head(res_Enz_DMSO$Symbol[which(res_Enz_DMSO$Group == 'Up-regulated')],15)
down.gene <- head(res_Enz_DMSO$Symbol[which(res_Enz_DMSO$Group == 'Down-regulated')],15)

res.top30.gene <- c(as.character(up.gene), as.character(down.gene))
res_Enz_DMSO$label[match(res.top30.gene, res_Enz_DMSO$Symbol)] <- res.top30.gene
tail(res_Enz_DMSO$logP,10)

ggscatter(res_Enz_DMSO, 
          x='log2FoldChange', y='logP',
          size = 1,
          label = res_Enz_DMSO$label,
          font.label = 8,
          color = "Group",
          palette = c('#0000FF','#BBBBBB','#FF0000'),
          xlab = 'log2FC',
          ylab = '-log10(p.adj)',
          repel = T) + theme_base()+xlim(-10,10)+ylim(0,400)


## ar sig
upgene <- res_Enz_DMSO[res_Enz_DMSO$Symbol[which(res_Enz_DMSO$log2FoldChange > 2)],]
upgene <- upgene[order(upgene$log2FoldChange, decreasing = T),]
ar_sig <- head(upgene$Symbol,50)
ar_heatmap <- assay(dds[ar_sig,])[,1:9]
ar_heatmap

annotation_col <- data.frame(Treatment = c("Veh","Veh","Veh","Enz","Enz","Enz","Enz*NRG1","Enz*NRG1","Enz*NRG1"))
rownames(annotation_col) <- colnames(ar_heatmap)
pheatmap(ar_heatmap, scale = 'row', color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)), 
         legend_labels='Z-score', cluster_cols=F, 
         annotation_col = annotation_col,main='AR Signature',
         angle_col = 45,fontsize_row=8) 
annotation_col
