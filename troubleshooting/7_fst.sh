#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name fst_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o fst.%j.out
#SBATCH -A lp_svbelleghem

#Step 1: map 100 resequencing short reads to reference primary genome (Dudzele and Nieuwpoort)

#Step 2: Calculate genotypes from 100 bam files, 1 VCF file with VCFtools

#Step 3: Strip VCF file into calls file (Simon Martin)

#Step 4: Calculate Fst values from the valls file, determine blocks of highest Fst values for identifying inversions
