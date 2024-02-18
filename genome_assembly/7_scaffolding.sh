#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name scaffolding
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o scaffolding.%j.out
#SBATCH -A lp_svbelleghem

#this is to tell my slurm job where conda is, so that it can activate my environment
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate thesis 

#YaHs will take the contig sequences (.fa) and the HiC aligned to the contigs (.bam produced in step 4) and scaffold them 
yahs /scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/Pogonus_hifiasm.asm.hic.p_ctg.fa \
/scratch/leuven/357/vsc35707/dudzele_pogonus/map_contigs/deduplicated_files/Pogonus_chalceus_r.bam \
-q 100000 -l 30 -r 50000

#final scaffolds (.fa) are used to run the ARIMA pipeline again, the mapping will be using the scaffolds as the reference (instead of the contigs)

module load SAMtools/1.13-GCC-10.3.0

#finally, the newly re-named assembly is also indexed in order to run the ARIMA pipeline once more, this time mapping to this scaffolding assembly
samtools faidx yahs.out_scaffolds_final.fa && cut -f1,2 yahs.out_scaffolds_final.fa.fai > yahs.out_scaffolds_final.fa.genome

#Can also scaffold using SALSA2
#bamToBed -i Pogonus_chalceus_r.bam > Pogonus_chalceus_r.bed
#sort -k 4 Pogonus_chalceus_r.bed > tmp && mv tmp Pogonus_chalceus_r.bed
#python run_pipeline.py -a Pogonus_hifiasm.asm.hic.p_ctg.fa -l Pogonus_hifiasm.asm.hic.p_ctg.fa.fai -b Pogonus_chalceus_r.bed -e DNASE -o scaffolds_salsa -m yes -g Pogonus_hifiasm.asm.hic.p_ctg.gfa
