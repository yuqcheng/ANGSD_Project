#! /bin/bash

#SBATCH --partition=panda
#SBATCH --nodes=1
#SBATCH --job-name=bamqc
#SBATCH --mem=16G

export PATH=$PATH:/softlib/apps/EL7/BamQC/bin

for i in *bam; do
bamqc -o /athena/angsd/scratch/yuc4009/bamqc -f /athena/angsd/scratch/yuc4009/refdata/hg38.ncbiRefSeq.gtf -g /athena/angsd/scratch/yuc4009/refdata -t 8 $i;
done
