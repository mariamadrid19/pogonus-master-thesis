#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name repeat_masker
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=48:00:00
#SBATCH -o masker.%j.out
#SBATCH -A lp_svbelleghem

conda activate EDTA2

RepeatMasker -species "tribolium castaneum" -pa 6 sorted_prim_dud.fasta
