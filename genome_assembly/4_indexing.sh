#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name indexing
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH -o indexing.%j.out
#SBATCH -A lp_svbelleghem

cd pogonus/hifiasm/

module load SAMtools/1.13-GCC-10.3.0
module load BWA/0.7.17-GCC-10.3.0

samtools faidx Pogonus_hifiasm.asm.hic.p_ctg.fa && cut -f1,2 Pogonus_hifiasm.asm.hic.p_ctg.fa.fai > Pogonus_hifiasm.asm.hic.p_ctg.fa.genome && bwa index Pogonus_hifiasm.asm.hic.p_ctg.fa
# this indexing is then used in the ARIMA mapping pipeline to generate the BAM file
