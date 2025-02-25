#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=20:00:00 
#SBATCH -A lp_svbelleghem

#Script written by Steven Van Belleghem (2024)

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BWA
module load SAMtools

echo "================="

# Sample IDs 
#samples=(\
#GC129388 GC129389 GC129390 GC129391 GC129392 \
#GC129393 GC129394 GC129395 GC129396 GC129397 \
#GC129398 GC129399 GC129400 GC129401 GC129402 \
#GC129403 GC129404 GC129405 GC129406 GC129407 \
#GC129408 GC129409 GC129410 GC129411 GC129412 \
#GC129413 GC129414 GC129415 GC129416 GC129417 \
#GC129418 GC129419 GC129420 GC129421 GC129422 \
#GC129423 GC129424 GC129425 GC129426 GC129427 \
#GC129428 GC129429 GC129430 GC129431 GC129432 \
#GC129433 GC129434 GC129435 GC129437 GC129438 \
#GC129439 GC129440 GC136078 GC136079 GC136080 \
#GC136081 GC136082 GC136083 GC136084 GC136085 \
#GC136086 GC136087 GC136088 GC136089 GC136090 \
#GC136091 GC136092 GC136093 GC136094 GC136095 \
#GC136096 GC136097 GC136098 GC136099 GC136100 \
#GC136101 GC136102 GC136103 GC136104 GC136105 \
#GC136106 GC136107 GC136108 GC136109 GC136110 \
#GC136111 GC136112 GC136113 GC136114 GC136115 \
#GC136116 GC136117 GC136118 GC136119 GC136120 \
#GC136121 GC136122 GC136123 GC136124 GC136125 \
#GC136126 GC136127 GC136128)

samples=(GC136096 GC136099)

echo "${samples[ID]}"

# Some folder and file paths to use later
#REF=/lustre1/scratch/350/vsc35085/Maria/sorted_prim_dud.fasta
REF=/lustre1/scratch/350/vsc35085/Maria/sorted_prim_nieu.fasta
#REFNAME=dudPrim
REFNAME=nieuPrim
BWAout=/lustre1/scratch/350/vsc35085/Maria/BAM

FILE1=/lustre1/scratch/350/vsc35085/Maria/Pogonus_RESEQ_ALL/$(echo "${samples[ID]}")_R1.fastq.gz
FILE2=/lustre1/scratch/350/vsc35085/Maria/Pogonus_RESEQ_ALL/$(echo "${samples[ID]}")_R2.fastq.gz

#if [ ! -f $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam ]; then

# Run BWA mapping
bwa mem -t 20 -M $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[ID]}").$REFNAME.bam

# Filter using samtools
samtools view -f 0x02 -q 20 -b $BWAout/$(echo "${samples[ID]}").$REFNAME.bam > $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam

# Sort using samtools
samtools sort $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam -o $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

# Remove PCR duplicates
java -jar /vsc-hard-mounts/leuven-data/350/vsc35085/programs/picard.jar MarkDuplicates \
-I $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam \
-O $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam \
-REMOVE_DUPLICATES true \
-M $BWAout/$(echo "${samples[ID]}").$REFNAME.dup_metrics.txt \
-ASSUME_SORTED true

# Remove intermediate files
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam 

samtools index $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam

cd /vsc-hard-mounts/leuven-data/350/vsc35085/scripts
sbatch -a 1-10 Run_mpileup.sh
#fi
echo "================="
