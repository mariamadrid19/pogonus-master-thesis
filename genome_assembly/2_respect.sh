#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name respect
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=24:00:00
#SBATCH -o respect.%j.out
#SBATCH -A lp_svbelleghem

conda activate respect

#tool to estimate genome size based on HiFi reads (fastq file)
respect -i /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/POG_HiFi_reads.fastq --threads 36
