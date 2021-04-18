#! /bin/bash

#SBATCH --partition=panda
#SBATCH --nodes=1
#SBATCH --job-name=bam_index
#SBATCH --mem=16G

spack load /qr4zqdd samtools@1.9%gcc@6.3.0 arch=linux-centos7-x86_64

for i in *.bam;do

echo "Start indexing: " $i;
samtools index $i;
done
