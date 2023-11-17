

juicer pre yahs.out.bin yahs.out_scaffolds_final.agp Ogibo_hifiasm2_4hap.asm.hic.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt
