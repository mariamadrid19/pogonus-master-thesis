#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name assemble_flye
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o hiflye.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

flye --pacbio-hifi /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/POG_HiFi_reads.fastq --genome-size 1.3g --threads 32
