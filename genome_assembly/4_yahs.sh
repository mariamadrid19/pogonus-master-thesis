#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name yahs
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=36:00:00
#SBATCH -o yahs.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/pogonus/yahs

module load SAMtools/1.13-GCC-10.3.0
module load BWA/0.7.17-GCC-10.3.0

samtools faidx /scratch/leuven/357/vsc35707/pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa && cut -f1,2 /scratch/leuven/357/vsc35707/pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa.fai > Pogonus_hifiasm.asm.hic.p_ctg.fa.genome && bwa index Pogonus_hifiasm.asm.hic.p_ctg.fa

yahs Pogonus_hifiasm.asm.hic.p_ctg.fa /scratch/leuven/357/vsc35707/pogonus/mapping/paired_bams/Pogonus_chalceus.bam -o results/POG
