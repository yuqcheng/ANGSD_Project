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

res_NRG1_DMSO <- results(dds, contrast = c("condition","DMSO", "NRG1"))
res_NRG1_DMSO <- as.data.frame(res_NRG1_DMSO)

res_NRG1_DMSO$Symbol <- rownames(res_NRG1_DMSO)

res_NRG1_DMSO <- res_NRG1_DMSO %>% dplyr::filter(abs(log2FoldChange)>0 & padj < 0.05)

res_NRG1_DMSO$logP <- -log10(res_NRG1_DMSO$padj)

res_NRG1_DMSO$Group <- 'not-significant'

res_NRG1_DMSO$Group[which(res_NRG1_DMSO$log2FoldChange > 2)] <- 'Down-regulated'
res_NRG1_DMSO$Group[which(res_NRG1_DMSO$log2FoldChange < -2)] <- 'Up-regulated'

res_NRG1_DMSO$label <- ""
res_NRG1_DMSO <- res_NRG1_DMSO[order(res_NRG1_DMSO$padj),]
up.gene <- head(res_NRG1_DMSO$Symbol[which(res_NRG1_DMSO$Group == 'Down-regulated')],15)
down.gene <- head(res_NRG1_DMSO$Symbol[which(res_NRG1_DMSO$Group == 'Up-regulated')],15)

res.top30.gene <- c(as.character(up.gene), as.character(down.gene))
res_NRG1_DMSO$label[match(res.top30.gene, res_NRG1_DMSO$Symbol)] <- res.top30.gene

ggscatter(res_NRG1_DMSO, 
          x='log2FoldChange', y='logP',
          size = 1,
          label = res_NRG1_DMSO$label,
          font.label = 8,
          color = "Group",
          palette = c('#0000FF','#BBBBBB','#FF0000'),
          xlab = 'log2FC',
          ylab = '-log10P',
          repel = T) + theme_base()+xlim(-10,10)+ylim(0,400)
