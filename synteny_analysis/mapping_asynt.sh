#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name asynt_map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o asynt_map.%j.out
#SBATCH -A lp_svbelleghem

module load SAMtools/1.13-GCC-10.3.0
conda activate thesis

samtools faidx scaffold_A_dud.fa 
samtools faidx scaffold_B_dud.fa 
samtools faidx scaffold_C_dud.fa 
samtools faidx scaffold_D_dud.fa 
samtools faidx scaffold_E_dud.fa 
samtools faidx scaffold_F_dud.fa
samtools faidx scaffold_G_dud.fa 
samtools faidx scaffold_H_dud.fa 

samtools faidx scaffold_A_nieu.fa
samtools faidx scaffold_B_nieu.fa
samtools faidx scaffold_C_nieu.fa
samtools faidx scaffold_D_nieu.fa
samtools faidx scaffold_E_nieu.fa
samtools faidx scaffold_F_nieu.fa
samtools faidx scaffold_G_nieu.fa
samtools faidx scaffold_H_nieu.fa


minimap2 -x asm10 scaffold_A_dud.fa scaffold_A_nieu.fa | gzip > scaffold_A.paf.gz

minimap2 -x asm10 scaffold_B_dud.fa scaffold_A_nieu.fa | gzip > scaffold_B.paf.gz

minimap2 -x asm10 scaffold_C_dud.fa scaffold_A_nieu.fa | gzip > scaffold_C.paf.gz

minimap2 -x asm10 scaffold_D_dud.fa scaffold_A_nieu.fa | gzip > scaffold_D.paf.gz

minimap2 -x asm10 scaffold_E_dud.fa scaffold_A_nieu.fa | gzip > scaffold_E.paf.gz

minimap2 -x asm10 scaffold_F_dud.fa scaffold_A_nieu.fa | gzip > scaffold_F.paf.gz

minimap2 -x asm10 scaffold_G_dud.fa scaffold_A_nieu.fa | gzip > scaffold_G.paf.gz

minimap2 -x asm10 scaffold_H_dud.fa scaffold_A_nieu.fa | gzip > scaffold_H.paf.gz
