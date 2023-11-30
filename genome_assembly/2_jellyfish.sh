#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name jellyfish
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH -o jellyfish.%j.out
#SBATCH -A lp_svbelleghem

#this will go into the seqtk directory to convert the fastq.gz reads into fasta reads
cd seqtk/
seqtk seq -a /scratch/leuven/357/vsc35707/pogonus/fastqc/POG_HiFi_reads.fastq.gz > /scratch/leuven/357/vsc35707/pogonus/fastqc/POG_HiFi_reads.fasta

module load Jellyfish/2.2.10-intel-2018a

#this will count kmers on the fasta file, takes a long time
jellyfish count -m 21 -s 100M -t 32 -C POG_HiFi_reads.fasta

#this makes the k-mer histogram (can be seen on the .out file, DO NOT DELETE IT)
jellyfish histo mer_counts.jf

#this generates an incredibly big file! be careful with storage
#jellyfish dump mer_counts.jf > mer_counts_dumps.fa

jellyfish info mer_counts.jf
