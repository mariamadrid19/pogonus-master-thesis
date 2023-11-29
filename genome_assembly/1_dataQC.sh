#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name fastqc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH -o fastqc.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

bam2fastq -o POG_HiFi_reads m64279e_231107_135307.reads.bam

module load FastQC/0.11.8-Java-1.8.0_162

fastqc POG_HiFi_reads.fastq.gz
