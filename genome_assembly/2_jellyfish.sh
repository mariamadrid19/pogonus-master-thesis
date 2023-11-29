#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name jellyfish
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=12:00:00
#SBATCH -o jellyfish.%j.out
#SBATCH -A lp_svbelleghem


module load Jellyfish/2.2.10-intel-2018a

jellyfish count -m 21 -s 100M -t 20 -C POG_HiFi_reads.fastq.gz

jellyfish histo mer_counts.jf

jellyfish dump mer_counts.jf > mer_counts_dumps.fa

jellyfish info mer_counts.jf
