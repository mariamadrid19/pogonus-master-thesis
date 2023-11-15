#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name rnaQUAST
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH -o rnaQUAST.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis #this is where rnaQUAST is installed
#https://github.com/ablab/rnaquast

rnaQUAST.py --transcripts /PATH/TO/transcripts1.fasta /PATH/TO/ANOTHER/transcripts2.fasta \
--busco insecta_lineage \
--gene_mark \
--o rnaQUAST_output
