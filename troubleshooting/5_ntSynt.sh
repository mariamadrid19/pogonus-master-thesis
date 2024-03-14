#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name synteny_ntsynt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -o ntsynt.%j.out
#SBATCH -A lp_svbelleghem

ntSynt dh1_50_scaffolds.fasta sorted_nh1.fasta -p ragtag_dh1_nh1 -d 20

python sort_ntsynt_blocks.py --synteny_blocks ragtag_dh1_nh1.synteny_blocks.tsv --sort_order dh1_50_scaffolds.fasta.fai sorted_nh1.fasta.fai --fais

python format_blocks_gggenomes.py --fai dh1_50_scaffolds.fasta.fai sorted_nh1.fasta.fai --prefix ragtag_dh1_nh1 --blocks ragtag_dh1_nh1.synteny_blocks.sorted.tsv --length 100 --colour dh1_50_scaffolds.fasta
