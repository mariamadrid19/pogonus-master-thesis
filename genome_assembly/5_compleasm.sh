#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name compleasm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH -o compleasm.%j.out
#SBATCH -A lp_svbelleghem

python compleasm.py run -a Pogonus_hifiasm.asm.hic.p_ctg.fa -o results/ -l insecta_db10 -t 32
