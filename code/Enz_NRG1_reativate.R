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
#DESeq.ds <- featurecount_ds

dds <- DESeq(DESeq.ds)

resultsNames(dds)

reactivate_gene <- results(dds, contrast = c("condition","DMSO", "Enz.NRG1"))
reactivate_gene <- subset(reactivate_gene, padj < 0.05 & abs(log2FoldChange) > 0)

rg <- as.data.frame(reactivate_gene)

rg$Symbol <- rownames(rg)

rg$logP <- -log10(rg$padj)

rg$Group <- 'not-significant'

rg$Group[which(rg$log2FoldChange > 2)] <- 'Down-regulated'
rg$Group[which(rg$log2FoldChange < -2)] <- 'Up-regulated'

rg$label <- ""
rg<- rg[order(rg$padj),]
up.gene <- head(rg$Symbol[which(rg$Group == 'Down-regulated')],15)
down.gene <- head(rg$Symbol[which(rg$Group == 'Up-regulated')],15)

res.top30.gene <- c(as.character(up.gene), as.character(down.gene))
rg$label[match(res.top30.gene, rg$Symbol)] <- res.top30.gene

rg$logP[which(rg$logP==Inf)] = 400

ggscatter(rg, 
          x='log2FoldChange', y='logP',
          size = 1,
          label = rg$label,
          font.label = 8,
          color = "Group",
          palette = c('#0000FF','#BBBBBB','#FF0000'),
          xlab = 'log2FC',
          ylab = '-log10P',
          repel = T) + 
  theme_base()+xlim(-10,10)+ylim(0,400) 

down <- length(which(reactivate_gene$log2FoldChange >0))
up <- length(which(reactivate_gene$log2FoldChange <0))

#reactivation sig

upgene <- rg[rg$Symbol[which(rg$log2FoldChange < -2)],]
upgene <- upgene[order(upgene$log2FoldChange, decreasing = F),]
re_sig <- head(upgene$Symbol,50)
re_sig <- c(re_sig, ar_sig)
re_heatmap <- assay(dds[rownames(rg),])
re_heatmap

annotation_col <- data.frame(Treatment = c("Veh","Veh","Veh","Enz","Enz","Enz","Enz*NRG1","Enz*NRG1","Enz*NRG1"))
rownames(annotation_col) <- colnames(re_heatmap[,1:9])
pheatmap(re_heatmap[,1:9], scale = 'row', color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)), 
         legend_labels='Z-score', cluster_cols=F, 
         annotation_col = annotation_col,main='Reactivation profile',
         angle_col = 45,show_rownames = F) 
