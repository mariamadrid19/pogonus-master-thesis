#! /bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name mapping_pipeline_juicer
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH -o juicer.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

cd /scratch/leuven/357/vsc35707/pogonus/yahs

(juicer pre yahs.out.bin yahs.out_scaffolds_final.agp Pogonus_hifiasm.asm.hic.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

java -jar -Xmx32G juicer_tools.1.9.9_jcuda.0.8.jar pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes) && (mv out.hic.part out.hic)
