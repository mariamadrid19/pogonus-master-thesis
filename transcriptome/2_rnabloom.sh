#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name rnabloom
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o rnabloom.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

rnabloom -long /scratch/leuven/357/vsc35707/pogonus/reads/Pogonus_PACBIO_RNA/POG_IsoSeq_HiFi_demux.fastq -lrpb -t 48 -outdir /scratch/leuven/357/vsc35707/pogonus/reads/Pogonus_PACBIO_RNA/TRANSCRIPTOME_ASSEMBLY/rnabloom_assembly
