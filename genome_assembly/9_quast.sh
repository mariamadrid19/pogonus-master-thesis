#! /bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name quast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH -o quast.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

quast /scratch/leuven/357/vsc35707/pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa -t 12