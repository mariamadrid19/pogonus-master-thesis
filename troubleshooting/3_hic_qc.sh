#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name hic_qc
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o indexing.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/

module load SAMtools/1.13-GCC-10.3.0
conda activate thesis #(bwa is loaded here)

bwa mem -5SP Pogonus_hifiasm.asm.hic.p_ctg.fa GC143248_ACTCTCGA-TGGTACAG_S65_R1.fastq GC143248_ACTCTCGA-TGGTACAG_S65_R2.fastq | samtools view -S -h -b -F 2304 > aligned_assembly_hic.bam

conda activate hic_qc

python hic_qc -b aligned_assembly_hic.bam -n 1000000 -r 
