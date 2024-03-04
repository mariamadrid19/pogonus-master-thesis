#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name assemble_canu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o hicanu.%j.out
#SBATCH -A lp_svbelleghem

canu -p hifi_canu_pogonus useGrid=false -d canu_results/ genomeSize=1300M -pacbio-hifi /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/POG_HiFi_reads.fastq
