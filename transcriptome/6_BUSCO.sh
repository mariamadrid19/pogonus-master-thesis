#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name BUSCO
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH -o busco.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis #this is where MetaEuk is installed)
module load BUSCO/3.0.2-foss-2018a-Python-2.7.14
module load HMMER/3.2.1-foss-2018a

busco --download-mode fast_lineage -l insecta -o insecta_lineage

busco -i assembled-transcriptome.fasta -l insecta -o busco_output  -m transcriptome
