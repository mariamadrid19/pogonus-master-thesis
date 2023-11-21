#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name rnaQUAST
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o rnaQUAST.%j.out
#SBATCH -A lp_svbelleghem

#this is where rnaQUAST is installed
conda activate thesis
#https://github.com/ablab/rnaquast

#https://cyverse.atlassian.net/wiki/spaces/TUT/pages/258736197/rnaQUAST+1.2.0+denovo+based+using+DE

rnaQUAST.py --transcripts /scratch/pogonus/reads/Pogonus_PACBIO_RNA/TRANSCRIPTOME_ASSEMBLY/BUSCO/ \
--busco insecta_lineage \
--o rnaQUAST_output_BUSCO 

rnaQUAST.py --transcripts /scratch/pogonus/reads/Pogonus_PACBIO_RNA/TRANSCRIPTOME_ASSEMBLY/BUSCO/ \
--gene_mark \
--o rnaQUAST_output_GM

#rnaQUAST.py --transcripts /PATH/TO/rnabloom_assembly.fasta /PATH/TO/isoscon_assembly.fasta /PATH/TO/isoseq3_assembly.fasta \
#--busco insecta_lineage \
#--o rnaQUAST_output_BUSCO

#rnaQUAST.py --transcripts /PATH/TO/rnabloom_assembly.fasta /PATH/TO/isoscon_assembly.fasta /PATH/TO/isoseq3_assembly.fasta \
#--gene_mark \
#--o rnaQUAST_output_GM

