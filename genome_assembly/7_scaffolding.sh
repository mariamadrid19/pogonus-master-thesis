#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name scaffolding
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o scaffolding.%j.out
#SBATCH -A lp_svbelleghem

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate thesis 

#YaHs will take the contig sequences (.fa) and the HiC aligned to the contigs (.bam produced in step 4) and scaffold them 
yahs /scratch/leuven/357/vsc35707/pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa /scratch/leuven/357/vsc35707/pogonus/map_contigs/deduplicated_files/Pogonus_chalceus_r.bam -o Pog_2.0

#Pog_2.0.fa is used to run the ARIMA pipeline again, the mapping will be using the scaffolds as the reference (instead of the contigs)

mv Pog_2.0_scaffolds_final.fa Pog_2.0.fa

module load SAMtools/1.13-GCC-10.3.0

#finally, the newly re-named assembly is also indexed in order to run the ARIMA pipeline once more, this time mapping to this scaffolding assembly
samtools faidx Pog_2.0.fa && cut -f1,2 Pog_2.0.fa.fai > Pog_2.0.fa.genome

