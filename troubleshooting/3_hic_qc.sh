#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name hic_qc
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o indexing.%j.out
#SBATCH -A lp_svbelleghem

module load SAMtools/1.13-GCC-10.3.0
module load BWA/0.7.17-GCC-10.3.0

bwa mem -5SP no_HiC_hifiasm.asm.bp.p_ctg.fa /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/GC143248_ACTCTCGA-TGGTACAG_S65_R1.fastq /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/GC143248_ACTCTCGA-TGGTACAG_S65_R2.fastq | samtools view -S -h -b -F 2304 > aligned_assembly_hic.bam

conda activate hic_qc

python hic_qc -b aligned_assembly_hic.bam -n 1000000 -r 
