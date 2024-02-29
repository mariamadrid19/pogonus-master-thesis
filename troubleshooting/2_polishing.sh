#Set input and parameters
round=2
threads=20
read=/scratch/leuven/357/vsc35707/dudzele_pogonus/hifiasm/POG_HiFi_reads.fastq
read_type=hifi
mapping_option=asm20
input=/scratch/leuven/357/vsc35707/dudzele_pogonus/HAP-1/purged_contigs_hap1.fasta

for ((i=1; i<=${round};i++)); do
    minimap2 -a -x ${mapping_option} -t ${threads} ${input} ${read}|samtools sort - -m 2g --threads ${threads} -o lgs.sort.bam;
    samtools index lgs.sort.bam;
    ls `pwd`/lgs.sort.bam > lgs.sort.bam.fofn;
    python NextPolish/lib/nextpolish2.py -g ${input} -l lgs.sort.bam.fofn -r ${read_type} -p ${threads} -sp -o genome.nextpolish.fa;
    if ((i!=${round}));then
        mv genome.nextpolish.fa genome.nextpolishtmp.fa;
        input=genome.nextpolishtmp.fa;
    fi;
done;
# Finally polished genome file: genome.nextpolish.fa
