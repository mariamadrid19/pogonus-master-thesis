#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name mpil 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10
#SBATCH --time=48:00:00 
#SBATCH -A lp_svbelleghem

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index  starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))


# Load the programs we will use
module load BWA
module load SAMtools
module load BCFtools
module load Python/3.7.0-foss-2018a

echo "================="

chrom=(CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag)
names=(1 10 2 3 4 5 6 7 8 9)


# Sample IDs 
samples=(\
GC129388 GC129389 GC129390 GC129391 GC129392 \
GC129393 GC129394 GC129395 GC129396 GC129397 \
GC129398 GC129399 GC129400 GC129401 GC129402 \
GC129403 GC129404 GC129405 GC129406 GC129407 \
GC129408 GC129409 GC129410 GC129411 GC129412 \
GC129413 GC129414 GC129415 GC129416 GC129417 \
GC129418 GC129419 GC129420 GC129421 GC129422 \
GC129423 GC129424 GC129425 GC129426 GC129427 \
GC129428 GC129429 GC129430 GC129431 GC129432 \
GC129433 GC129434 GC129435 GC129437 GC129438 \
GC129439 GC129440 GC136078 GC136079 GC136080 \
GC136081 GC136082 GC136083 GC136084 GC136085 \
GC136086 GC136087 GC136088 GC136089 GC136090 \
GC136091 GC136092 GC136093 GC136094 GC136095 \
GC136096 GC136097 GC136098 GC136099 GC136100 \
GC136101 GC136102 GC136103 GC136104 GC136105 \
GC136106 GC136107 GC136108 GC136109 GC136110 \
GC136111 GC136112 GC136113 GC136114 GC136115 \
GC136116 GC136117 GC136118 GC136119 GC136120 \
GC136121 GC136122 GC136123 GC136124 GC136125 \
GC136126 GC136127 GC136128)

# Some folder and file paths to use later
#REF=/lustre1/scratch/350/vsc35085/Maria/sorted_prim_dud.fasta
REF=/lustre1/scratch/350/vsc35085/Maria/sorted_prim_nieu.fasta
#REFNAME=dudPrim
REFNAME=nieuPrim
BWAout=/lustre1/scratch/350/vsc35085/Maria/BAM

# make a single list of all the samples that can be used in the samtools command
ALL_LIST=""
for FILE in ${samples[*]}
do
ALL_LIST="$ALL_LIST $FILE".$REFNAME.filtered.sorted.nd.bam""
done
eval command=\$$(echo ALL_LIST)

# run mpileup
cd /lustre1/scratch/350/vsc35085/Maria/BAM

#bcftools mpileup -O z --threads 10 -f $REF $(echo $command) -r $(echo "${chrom[ID]}") | bcftools call -m -Oz -o /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").vcf.gz 

module load VCFtools
module load tabix

vcftools --gzvcf /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").vcf.gz --recode --remove-indels --minQ 30 --max-missing 0.25 --stdout | bgzip > /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").filt.vcf.gz

#bgzip /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").filt.vcf

python /vsc-hard-mounts/leuven-data/350/vsc35085/programs/parseVCF.py \
--minQual 30 \
--skipIndels \
-i /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").filt.vcf.gz  | \
gzip > /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").calls.gz

zcat /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").calls.gz | sed 's/.nieuPrim.filtered.sorted.nd.bam//g' | bgzip > /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").H.calls.gz

#Script written by Steven Van Belleghem (2024)
