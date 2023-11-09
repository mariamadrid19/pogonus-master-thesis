#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name rnabloom
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o rnabloom.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

rnabloom -long POG_larveIsoSeq_HiFi.fastq -t 30 -lrpb -outdir /scratch/leuven/357/vsc35707/pogonus/reads/Pogonus_PACBIO_RNA/rnabloom_assembly
