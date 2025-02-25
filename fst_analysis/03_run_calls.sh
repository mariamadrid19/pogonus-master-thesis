#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name calls 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00 
#SBATCH -A lp_svbelleghem

#Script written by Steven Van Belleghem (2024)


ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load Python/3.7.0-foss-2018a
#export PYTHONPATH=$PYTHONPATH:/vsc-hard-mounts/leuven-data/350/vsc35085/programs/genomics_general
module load tabix

#module load BCFtools

echo "================="

chrom=(CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag)
names=(1 10 2 3 4 5 6 7 8 9)
#REFNAME=dudPrim
REFNAME=nieuPrim

#bcftools view /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.fix.vcf.gz --regions $(echo "${chrom[ID]}") | gzip > /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").vcf.gz

vcftools --gzvcf /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").vcf.gz --recode --remove-indels --minQ 30 --max-missing 0.25 --stdout | bgzip > /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").filt.vcf.gz

python /vsc-hard-mounts/leuven-data/350/vsc35085/programs/parseVCF.py \
--minQual 30 \
--skipIndels \
-i /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").filt.vcf.gz  | \
gzip > /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").calls.gz


zcat /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").calls.gz | sed 's/.nieuPrim.filtered.sorted.nd.bam//g' | gzip > /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").H.calls.gz


#--gtf flag=GQ min=30 gtTypes=Het \
#--gtf flag=GQ min=30 gtTypes=HomAlt \
#--gtf flag=DP min=10 max=100 \
