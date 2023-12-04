#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name gfastats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH -o gfastats.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

gfastats /scratch/leuven/357/vsc35707/pogonus/map_scaffolds/yahs.out_scaffolds_final.fa --threads 32
