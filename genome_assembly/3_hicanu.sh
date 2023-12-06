#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name assemble_canu
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o hicanu.%j.out
#SBATCH -A lp_svbelleghem

export PATH=$PATH:/scratch/leuven/357/vsc35707/pogonus/canu-2.2/build/bin/

canu -p hifi_canu_pogonus -d canu_results/ genomeSize=444m -pacbio-hifi /scratch/leuven/357/vsc35707/pogonus/hifiasm/POG_HiFi_reads.fastq.gz 
