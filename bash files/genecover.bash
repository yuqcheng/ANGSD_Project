#! /bin/bash

#SBATCH --partition=panda
#SBATCH --job-name=rseqc
#SBATCH --nodes=1
#SBATCH --mem=16G

echo "Start RSeQC"
for i in /athena/angsd/scratch/yuc4009/bams_selected/*.bam; do
echo '$i';
geneBody_coverage.py -r /athena/angsd/scratch/yuc4009/refdata/hg38_RefSeq.bed -i $i -o $i;
done

spack load -r py-multiqc@1.7

multiqc .
