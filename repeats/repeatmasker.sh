#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name repeat_masker
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=48:00:00
#SBATCH -o masker.%j.out
#SBATCH -A lp_svbelleghem

conda activate repeats
 
BuildDatabase -name dudRM sorted_prim_dud.fasta
RepeatModeler -database dudRM -threads 24 -LTRStruct >& run.out
RepeatMasker -lib dudRM-families.fa -pa 24 sorted_prim_dud.fasta -dir repeat_trial2 -gff
 
BuildDatabase -name nieRM sorted_prim_nieu.fasta
RepeatModeler -database nieRM -threads 24 -LTRStruct >& run2.out
RepeatMasker -lib nieM-families.fa -pa 24 sorted_prim_nieu.fasta -dir repeat_trial_nie -gff
