#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name isoseq3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o isoseq3.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

#clustering step
isoseq3 cluster POG_larveIsoSeq.demux.bam POG_larveIsoSeq.unpolished.bam --verbose
#polishing step is not necessary (subreads are not available)

isoseq3 summarize POG_larveIsoSeq.unpolished.bam summary-isoseq3.csv

