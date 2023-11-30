#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name jellyfish
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH -o jellyfish.%j.out
#SBATCH -A lp_svbelleghem

cd seqtk/
seqtk seq -a /scratch/leuven/357/vsc35707/pogonus/fastqc/POG_HiFi_reads.fastq.gz > /scratch/leuven/357/vsc35707/pogonus/fastqc/POG_HiFi_reads.fasta

module load Jellyfish/2.2.10-intel-2018a

jellyfish count -m 21 -s 100M -t 32 -C POG_HiFi_reads.fasta

jellyfish histo mer_counts.jf

jellyfish dump mer_counts.jf > mer_counts_dumps.fa

jellyfish info mer_counts.jf
