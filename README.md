# ANGSD_project
This is used for the final project of CMPB 5004 in Weill Cornell.

This project mainly focus on the gene expression analysis of CRPC.

Androgen receptor (AR) is a survival factor for prostate cancer cells that plays an essential role in cancer progression. However, after a period of androgen deprivation treatment (ADT), cancer cells will develop resistance to castration. Previous studies have shown that the acquisition of such resistance is caused by the tumor microenvironment. A factor called NRG1 plays an important role in the resistance gaining process. In this study, we hypothesize that castration-resistant cancer cells (CRPCs) induced by NRG1 have different gene expression profiles from normal prostate cancer cells. This drug resistance is not entirely caused by the reactivation of AR-related pathways. There might be other pathways that allow cancer cells to survive, and CRPCs have different activated functional clusters from prostate cancer cells that have not developed drug resistance.

Here, we show our analysis pipeline as follows.

![stat](https://github.com/yuqcheng/ANGSD_project/blob/main/figure/pipeline.png)

We collected the bulk RNA-seq data of cancer cells with different treatments. We choose the sample treated with vehicle as control, the sample treated with Enzalutamide (Enz) as AR positive experiment and the sample treated with Enz and NRG1 and the reactivation experiment. 

## Samples

RNA-Seq libraries were prepared using the Illumina TruSeq stranded mRNA kit, with 10 cycles of PCR amplification, starting from 500 ng of total RNA, at the Genome Technology Center (GTC) at New York University Langone Medical Center.

Here, we have 4 samples (Cancer cells with different treatment), each sample has 3 replicates. 
- **Sample 1: DMSO (treated with DMSO (vehicle))**
- **Sample 2: Enz (treated with Enzalutamide, a androgen receptor (AR) target inhibitor)**
- **Sample 3: NRG1 (treated with recombinant NRG1 peptide, a factor that induce CRPCs gain resistance to Enz)**
- **Sample 4: Enz+NRG1 (treated with both Enz and NRG1 to show the resistance to Enz)**

Here you can find all my fastq files and bam files in my scratch folder. (Please note that all the fastq files are pair-ended)

**Fastq PATH: /athena/angsd/scratch/yuc4009/NRG1**

**Trimmed fastq: /athena/angsd/scratch/yuc4009/NRG1/trim**

**Bams PATH: /athena/angsd/scratch/yuc4009/bams && /athena/angsd/scratch/yuc4009/bams_selected**

**Bamqc PATH: /athena/angsd/scratch/yuc4009/bamqc**

In this repo, we show all of the bash files that we used to run STAR.

## STAR

We construted the index folder with hg38.fa and hg38.ncbirefdata.gz. All the bash files are uploaded in the bashs file folder.

**STAR index PATH: /athena/angsd/scratch/yuc4009/refdata**

## FeatureCounts

We use all 4\*3=12 bam files (collected in the bams_selected) to run featurecounts. The result is shown in the **/athena/angsd/scratch/yuc4009/featureCounts**.

In this repo, we show the featureCounts result (.txt) and exploratory analysis result (.Rmd)

Here we show the stat of the featureCounts result.

![stat](https://github.com/yuqcheng/ANGSD_project/blob/main/figure/featurecounts_stat.png)

Then we generated dendrogram and PCA plots. The dendrogram is shown in the Rmd, we only show the PCA plot here. We can clearly see that the pca result is correspond with our 4 samples.

![pca](https://github.com/yuqcheng/ANGSD_project/blob/main/figure/pca.png)

From the PCA plot, we can easily find that the 4 samples are clustered well, which means these four samples have a quite large difference inter-sample while a very small difference intra-replicates. This conclusion is biological meaningfully. 
- For DMSO, it is only treated with the vehicle so that it should have a normal AR activity. 
- For Enz, the sample is treated with Enzalutamide so that its AR activity is inhibited and the expression of AR-target genes are changes. 
- For NRG1, the sample is treated with the recombinant NRG1 so that the NRG1-HER2 pathway is activated and NRG1-HER2 targeted genes are activated/inhibited. 
- For Enz+NRG1, the AR pathway is inhibited and NRG1-HER2 pathway is activated. In this case, we can find what happens after the cancer cells gain drug resistance.

So it is natural that we can gain 4 distinct clusters after PCA.
