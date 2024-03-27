#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name hap_BWA_map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH -o hap_BWA_map.%j.out
#SBATCH -A lp_svbelleghem

module load BWA/0.7.17-GCC-10.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.18.23-Java-1.8.0_171

indID=$((SLURM_ARRAY_TASK_ID -1))

# Sample IDs
samples=(\
GC129388 GC129389 \
GC129394 GC129395)

echo "${samples[indID]}"

REF=/scratch/leuven/357/vsc35707/dudzele_pogonus/HAP-1/yahs.out_scaffolds_final.fa
READS=/scratch/leuven/357/vsc35707/Pogonus_RESEQ_ALL/
BWAout=/scratch/leuven/357/vsc35707/Pogonus_RESEQ_ALL/bams_coverage_analysis_dh1

FILE1=$READS/$(echo "${samples[indID]}")
FILE2=$READS/$(echo "${samples[indID]}")

# get the complete filename
FILE1="$(ls  $READS | grep "$(echo "${samples[indID]}")" | grep "R1")"
FILE2="$(ls  $READS | grep "$(echo "${samples[indID]}")" | grep "R2")"

# Run BWA mapping
echo "### Step 1: FASTQ to BAM"
bwa mem -t 8 -k 17 -M $REF $READS/$FILE1 $READS/$FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[indID]}").Pchal_DH1.bam

# Filter using samtools
echo "### Step 2: FILTER using samtools"
samtools view -f 0x02 -q 20 -b $BWAout/GC129388.Pchal_DH1.bam > $BWAout/GC129388.Pchal_DH1.filtered.bam

# Sort using samtools
echo "### Step 3: SORT using samtools"
samtools sort $BWAout/$(echo "${samples[indID]}").Pchal_DH1.filtered.bam -o $BWAout/$(echo "${samples[indID]}").Pchal_DH1.filtered.sorted.bam

# Remove PCR duplicates
echo "### Step 4: REMOVE PCR duplicates using PICARD"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$BWAout/$(echo "${samples[indID]}").Pchal_DH1.filtered.sorted.bam \
OUTPUT=$BWAout/$(echo "${samples[indID]}").Pchal_DH1.filtered.sorted.nd.bam REMOVE_DUPLICATES=true \
METRICS_FILE=$BWAout/$(echo "${samples[indID]}").Pchal_DH1.dup_metrics.txt ASSUME_SORTED=true

# Remove intermediate files
echo "### Step 5: REMOVE intermediate files"
#rm $BWAout/$(echo "${samples[indID]}").Pchal_DH1.bam
rm $BWAout/$(echo "${samples[indID]}").Pchal_DH1.filtered.bam
rm $BWAout/$(echo "${samples[indID]}").Pchal_DH1.filtered.sorted.bam
