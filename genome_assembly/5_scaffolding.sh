#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name scaffolding
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=36:00:00
#SBATCH -o scaffolding.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/pogonus/yahs

module load SAMtools/1.13-GCC-10.3.0
module load BWA/0.7.17-GCC-10.3.0

samtools faidx /scratch/leuven/357/vsc35707/pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa && cut -f1,2 /scratch/leuven/357/vsc35707/pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa.fai > Pogonus_hifiasm.asm.hic.p_ctg.fa.genome && bwa index Pogonus_hifiasm.asm.hic.p_ctg.fa

yahs Pogonus_hifiasm.asm.hic.p_ctg.fa /scratch/leuven/357/vsc35707/pogonus/mapping/deduplicated_files/Pogonus_chalceus_rep1.bam -o results/POG

#the yahs.out_scaffolds_final.fa file produced here (in the POG directory) is then re-named Pog_2.0.fasta and used to run the ARIMA pipeline again
#this time, the mapping will be using the Pog_2.0.fasta as the reference (instead of the hifiasm assembly)

cp POG/yahs.out_scaffolds_final.fa /scratch/leuven/357/vsc35707/pogonus/mapping_HiC
cd /scratch/leuven/357/vsc35707/pogonus/mapping_HiC
mv yahs.out_scaffolds_final.fa Pog_2.0.fa

samtools faidx Pog_2.0.fa && cut -f1,2 Pog_2.0.fa.fai > Pog_2.0.fa.genome && bwa index Pog_2.0.fa

