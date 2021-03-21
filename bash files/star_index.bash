#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=star_index
#SBATCH --mem=50G

STAR --runMode genomeGenerate --runThreadN 8 --genomeDir star_index --genomeFastaFiles /athena/angsd/scratch/yuc4009/refdata/hg38.fa --sjdbGTFfile /athena/angsd/scratch/yuc4009/refdata/hg38.ncbiRefSeq.gtf --sjdbOverhang 99
