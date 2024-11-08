#Script written by Steven Van Belleghem (2024)

#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10 
#SBATCH --time=20:00:00 
#SBATCH -A lp_svbelleghem

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load Python/3.7.0-foss-2018a
module load tabix
export PYTHONPATH=$PYTHONPATH:/vsc-hard-mounts/leuven-data/350/vsc35085/programs/genomics_general

REFNAME=dudPrim
#REFNAME=nieuPrim

names=(1 10 2 3 4 5 6 7 8 9)

zcat /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").calls.gz | sed 's/.dudPrim.filtered.sorted.nd.bam//g' | bgzip > /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$(echo "${names[ID]}").H.calls.gz
