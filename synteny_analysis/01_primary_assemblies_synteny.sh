#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name contigous_scaffolds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o contigous_scaffolds.%j.out
#SBATCH -A lp_svbelleghem

conda activate ragtag

#First step is to sort the primary assemblies based on the old reference assembly (biggest 50 scaffolds)
ragtag.py scaffold 50_ref_scaffolds.fasta primary_dudzele.fa -t 24 -o dud_prim

ragtag.py scaffold 50_ref_scaffolds.fasta primary_nieuwpoort.fa -t 24 -o nieu_prim

#Then, the names of the RagTag files need to be changed to obtain the sorted and merged scaffolds for each of the primary assemblies

#With that, the 10 linkage groups are extracted from each assembly. This way, we obtain the 10 "chromosomes" for each of the assemblies
awk '/^>/ { if (count++ >= 10) exit } { print }' sorted_primary_dudzele.fasta > 10_sorted_prim_dud.fasta

awk '/^>/ { if (count++ >= 10) exit } { print }' sorted_primary_nieuwpoort.fasta > 10_sorted_prim_nieu.fasta

#We can then study the synteny between them 
conda activate ntsynt

ntSynt 10_sorted_prim_dud.fasta 10_sorted_prim_nieu.fasta -p primary_10 -t 24 -d 30

python denovo_synteny_block_stats.py --tsv primary_10.synteny_blocks.tsv --fai 10_sorted_prim_dud.fasta.fai 10_sorted_prim_nieu.fasta.fai

python sort_ntsynt_blocks.py --synteny_blocks primary_10.synteny_blocks.tsv --sort_order 10_sorted_prim_dud.fasta.fai 10_sorted_prim_nieu.fasta.fai --fais > primary_10.synteny_blocks.sorted.tsv

python format_blocks_gggenomes.py --fai 10_sorted_prim_dud.fasta.fai 10_sorted_prim_nieu.fasta.fai --prefix primary_10 --blocks primary_10.synteny_blocks.sorted.tsv --length 100 --colour 10_sorted_prim_dud.fasta

cp primary_10.links.tsv $VSC_DATA
cp primary_10.sequence_lengths.tsv $VSC_DATA
