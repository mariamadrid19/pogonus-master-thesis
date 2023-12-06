#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name minimap2
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o map_scaffolds_ref.%j.out
#SBATCH -A lp_svbelleghem

module load minimap2/2.12-foss-2018a

#https://www.biostars.org/p/9552255/

#index reference genome
minimap2 -d GCA_002278615.1_Pchal_1.0_genomic.mmi GCA_002278615.1_Pchal_1.0_genomic.fna

#map the scaffolds to the reference genome
minimap2 -ax asm5 -t 32 GCA_002278615.1_Pchal_1.0_genomic.fna 50_scaffolds.fa > mapped_scaffolds_ref.paf
