#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name purge_dups
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o purgedups.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

#Split an assembly and do a self-self alignment
split_fa Pogonus_hifiasm.asm.hic.p_ctg.fa > Pogonus_hifiasm.asm.hic.p_ctg.fa.split

minimap2 -xmap-hifi -DP Pogonus_hifiasm.asm.hic.p_ctg.fa.split Pogonus_hifiasm.asm.hic.p_ctg.fa.split | gzip -c - > Pogonus_hifiasm.asm.hic.p_ctg.fa.split.self.paf.gz

#purge haplotigs and overlaps 
purge_dups -2 -T cutoffs -c PB.base.cov Pogonus_hifiasm.asm.hic.p_ctg.fa.split.self.paf.gz > dups.bed 2> purge_dups.log

#get purged primary and haplotig sequences from draft assembly
get_seqs -e dups.bed Pogonus_hifiasm.asm.hic.p_ctg.fa
