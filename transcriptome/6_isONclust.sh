#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name isonclust
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o isonclust.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

#clustering step
isONclust --isoseq --fastq POG_IsoSeq_HiFi_demux.fastq --outfolder TRANSCRIPTOME_ASSEMBLY/isONclust/ --t 32
