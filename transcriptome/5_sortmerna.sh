#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name sortmerna
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o sortmerna.%j.out

conda activate thesis
cd /scratch/leuven/357/vsc35707/pogonus/reads/Pogonus_PACBIO_RNA/rRNA

sortmerna --ref silva_rfam_databases/rfam-5.8s-database-id98.fasta --ref silva_rfam_databases/rfam-5s-database-id98.fasta --ref silva_rfam_databases/silva-euk-18s-id95.fasta --ref silva_rfam_databases/silva-euk-28s-id98.fasta --reads POG_IsoSeq_HiFi_demux.fasta
