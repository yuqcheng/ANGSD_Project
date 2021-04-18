#! /bin/bash

#SBATCH --partition=panda
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --job-name=fastqc

source ~/.bash_profile^C
spack load fastqc

for i in {57..77}; do
fastqc -o /athena/angsd/scratch/yuc4009/fastqc -t 6 /athena/angsd/scratch/yuc4009/NRG1/SRR114706${i}_1.fastq.gz /athena/angsd/scratch/yuc4009/NRG1/SRR114706${i}_2.fastq.gz;
done;
