#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name dh2_BWA_map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o dh2_BWA_map.%j.out
#SBATCH -A lp_edu_comparativegenomics

module load BWA/0.7.17-GCC-10.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.18.23-Java-1.8.0_171

REF=/scratch/leuven/357/vsc35707/dudzele_pogonus/HAP-1/yahs.out_scaffolds_final.fa
READS=/scratch/leuven/357/vsc35707/Pogonus_RESEQ_ALL/
BWAout=/scratch/leuven/357/vsc35707/Pogonus_RESEQ_ALL/bams_coverage_analysis_dh1

echo "### READ PAIR 1"
# READS 1
# Run BWA mapping
echo "### Step 1: FASTQ to BAM"
bwa mem -t 8 -k 17 -M $REF $READS/GC129388_R1.fastq.gz $READS/GC129388_R2.fastq.gz | samtools view -bS - > $BWAout/GC1
29388.Pchal_DH1.bam

# Filter using samtools
echo "### Step 2: FILTER using samtools"
samtools view -f 0x02 -q 20 -b $BWAout/GC129388.Pchal_DH1.bam > $BWAout/GC129388.Pchal_DH1.filtered.bam

# Sort using samtools
echo "### Step 3: SORT using samtools"
samtools sort $BWAout/GC129388.Pchal_DH1.filtered.bam -o $BWAout/GC129388.Pchal_DH1.filtered.sorted.bam

# Remove PCR duplicates
echo "### Step 4: REMOVE PCR duplicates using PICARD"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$BWAout/GC129388.Pchal_DH1.filtered.sorted.bam OUTPUT=$BWAout/G
C129388.Pchal_DH1.filtered.sorted.nd.bam REMOVE_DUPLICATES=true METRICS_FILE=$BWAout/GC129388.Pchal_DH1.dup_metrics.txt
 ASSUME_SORTED=true

# Remove intermediate files
echo "### Step 5: REMOVE intermediate files"
#rm $BWAout/GC129388.Pchal_DH1.bam
#rm $BWAout/GC129388.Pchal_DH1.filtered.bam
#rm $BWAout/GC129388.Pchal_DH1.filtered.sorted.bam

echo "### READ PAIR 2"
# READS 2
# Run BWA mapping
echo "### Step 1: FASTQ to BAM"
bwa mem -t 8 -k 17 -M $REF $READS/GC129389_R1.fastq.gz $READS/GC129389_R2.fastq.gz | samtools view -bS - > $BWAout/GC1
29389.Pchal_DH1.bam

# Filter using samtools
echo "### Step 2: FILTER using samtools"
samtools view -f 0x02 -q 20 -b $BWAout/GC129389.Pchal_DH1.bam > $BWAout/GC129389.Pchal_DH1.filtered.bam

# Sort using samtools
echo "### Step 3: SORT using samtools"
samtools sort $BWAout/GC129389.Pchal_DH1.filtered.bam -o $BWAout/GC129389.Pchal_DH1.filtered.sorted.bam

# Remove PCR duplicates
echo "### Step 4: REMOVE PCR duplicates using PICARD"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$BWAout/GC129389.Pchal_DH1.filtered.sorted.bam OUTPUT=$BWAout/G
C129389.Pchal_DH1.filtered.sorted.nd.bam REMOVE_DUPLICATES=true METRICS_FILE=$BWAout/GC129389.Pchal_DH1.dup_metrics.txt
 ASSUME_SORTED=true

# Remove intermediate files
echo "### Step 5: REMOVE intermediate files"
#rm $BWAout/GC129389.Pchal_DH1.bam
#rm $BWAout/GC129389.Pchal_DH1.filtered.bam
#rm $BWAout/GC129389.Pchal_DH1.filtered.sorted.bam

echo "### READ PAIR 3"
# READS 3
# Run BWA mapping
echo "### Step 1: FASTQ to BAM"
bwa mem -t 8 -k 17 -M $REF $READS/GC129394_R1.fastq.gz $READS/GC129394_R2.fastq.gz | samtools view -bS - > $BWAout/GC1
29394.Pchal_DH1.bam

# Filter using samtools
echo "### Step 2: FILTER using samtools"
samtools view -f 0x02 -q 20 -b $BWAout/GC129394.Pchal_DH1.bam > $BWAout/GC129394.Pchal_DH1.filtered.bam

# Sort using samtools
echo "### Step 3: SORT using samtools"
samtools sort $BWAout/GC129394.Pchal_DH1.filtered.bam -o $BWAout/GC129394.Pchal_DH1.filtered.sorted.bam

# Remove PCR duplicates
echo "### Step 4: REMOVE PCR duplicates using PICARD"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$BWAout/GC129394.Pchal_DH1.filtered.sorted.bam OUTPUT=$BWAout/G
C129394.Pchal_DH1.filtered.sorted.nd.bam REMOVE_DUPLICATES=true METRICS_FILE=$BWAout/GC129394.Pchal_DH1.dup_metrics.txt
 ASSUME_SORTED=true

# Remove intermediate files
echo "### Step 5: REMOVE intermediate files"
#rm $BWAout/GC129394.Pchal_DH1.bam
#rm $BWAout/GC129394.Pchal_DH1.filtered.bam
#rm $BWAout/GC129394.Pchal_DH1.filtered.sorted.bam

echo "### READ PAIR 4"
# READS 4
# Run BWA mapping
echo "### Step 1: FASTQ to BAM"
bwa mem -t 8 -k 17 -M $REF $READS/GC129395_R1.fastq.gz $READS/GC129395_R2.fastq.gz | samtools view -bS - > $BWAout/GC12
9395.Pchal_DH1.bam

# Filter using samtools
echo "### Step 2: FILTER using samtools"
samtools view -f 0x02 -q 20 -b $BWAout/GC129395.Pchal_DH1.bam > $BWAout/GC129395.Pchal_DH1.filtered.bam

# Sort using samtools
echo "### Step 3: SORT using samtools"
samtools sort $BWAout/GC129395.Pchal_DH1.filtered.bam -o $BWAout/GC129395.Pchal_DH1.filtered.sorted.bam

# Remove PCR duplicates
echo "### Step 4: REMOVE PCR duplicates using PICARD"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$BWAout/GC129395.Pchal_DH1.filtered.sorted.bam OUTPUT=$BWAout/G
C129395.Pchal_DH1.filtered.sorted.nd.bam REMOVE_DUPLICATES=true METRICS_FILE=$BWAout/GC129395.Pchal_DH1.dup_metrics.txt
 ASSUME_SORTED=true

# Remove intermediate files
echo "### Step 5: REMOVE intermediate files"
#rm $BWAout/GC129395.Pchal_DH1.bam
#rm $BWAout/GC129395.Pchal_DH1.filtered.bam
#rm $BWAout/GC129395.Pchal_DH1.filtered.sorted.bam
