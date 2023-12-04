#! /bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name map_contigs
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=16
#SBATCH --time=48:00:00
#SBATCH -o map_contigs.%j.out
#SBATCH -A lp_svbelleghem

conda activate omniC

bwa mem -5SP -T0 -t 16 /scratch/leuven/357/vsc35707/pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa /scratch/leuven/357/vsc35707/pogonus/hifiasm/GC143248_ACTCTCGA-TGGTACAG_S65_R1.fastq /scratch/leuven/357/vsc35707/pogonus/hifiasm/GC143248_ACTCTCGA-TGGTACAG_S65_R2.fastq | \
pairtools parse --min-mapq 40 --walks-policy 5unique \
--max-inter-align-gap 30 --nproc-in 16 --nproc-out 16 --chroms-path /scratch/leuven/357/vsc35707/pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa.genome | \
pairtools sort --tmpdir=/temp/ --nproc 12|pairtools dedup --nproc-in 16 \
--nproc-out 16 --mark-dups --output-stats stats.txt|pairtools split --nproc-in 16 \
--nproc-out 16 --output-pairs mapped.pairs --output-sam -|samtools view -bS -@16 | \
samtools sort -@16 -o mapped.PT.bam;samtools index mapped.PT.bam
