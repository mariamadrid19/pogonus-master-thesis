#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name assemble_canu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o hicanu.%j.out
#SBATCH -A lp_svbelleghem
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=20000M

#Canu installation
#curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.tar.xz --output canu-2.2.tar.xz
#tar -xJf canu-2.2.tar.xz
#cd canu-2.2/src
#make -j 8

#export PATH=$PATH:/scratch/leuven/357/vsc35707/dudzele_pogonus/canu-2.2/build/bin/

canu -p hifi_canu_pogonus -d canu_results/ genomeSize=1300M -pacbio-hifi /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/POG_HiFi_reads.fastq 
