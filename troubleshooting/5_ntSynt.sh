#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name synteny_ntsynt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -o ntsynt_all.%j.out
#SBATCH -A lp_svbelleghem

conda activate ntsynt 

ntSynt dh1_50_scaffolds.fasta sorted_dh2.fasta sorted_nh1.fasta sorted_nh2.fasta -p synteny_all_haps -t 24 -d 30

python denovo_synteny_block_stats.py --tsv synteny_all_haps.synteny_blocks.tsv --fai dh1_50_scaffolds.fasta.fai sorted_dh2.fasta.fai \
sorted_nh1.fasta.fai sorted_nh2.fasta.fai

python sort_ntsynt_blocks.py --synteny_blocks synteny_all_haps.synteny_blocks.tsv --sort_order dh1_50_scaffolds.fasta.fai \
sorted_dh2.fasta.fai sorted_nh1.fasta.fai sorted_nh2.fasta.fai --fais > synteny_all_haps.synteny_blocks.sorted.tsv

python format_blocks_gggenomes.py --fai dh1_50_scaffolds.fasta.fai sorted_dh2.fasta.fai sorted_nh1.fasta.fai sorted_nh2.fasta.fai \
--prefix synteny_all_haps --blocks synteny_all_haps.synteny_blocks.sorted.tsv --length 100 --colour dh1_50_scaffolds.fasta

cat synteny_all_haps.links.tsv  | mlr --tsv sort -f strand -n block_id > synteny_all_haps.links.sorted.tsv && mv synteny_all_haps.links.sorted.tsv synteny_all_haps.links.tsv
cat synteny_all_haps.sequence_lengths.tsv | mlr --tsv sort -f seq_id > synteny_all_haps.sequence_lengths.sorted.tsv && mv synteny_all_haps.sequence_lengths.sorted.tsv synteny_all_haps.sequence_lengths.tsv

#run R script to generate plots
Rscript plot_synteny_blocks_gggenomes.R -s synteny_all_haps.sequence_lengths.tsv -l synteny_all_haps.links.tsv --scale 50000000 --p synteny_all_haps
