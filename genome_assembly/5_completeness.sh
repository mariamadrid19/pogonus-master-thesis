#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name compleasm_quast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH -o compleasm_quast.%j.out
#SBATCH -A lp_svbelleghem

# download compleasm and its dependencies (miniprot and hmmsearch)
#wget https://github.com/huangnengCSU/compleasm/releases/download/v0.2.5/compleasm-0.2.5_x64-linux.tar.bz2
#tar -jxvf compleasm-0.2.5_x64-linux.tar.bz2

#conda install pandas

# Run compleasm if lineage is known
compleasm_kit/compleasm.py download insecta
compleasm_kit/compleasm.py run -a /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa -o results/ -l insecta -t 32

conda activate thesis

quast /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa -t 12
