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
#module load tabix
export PYTHONPATH=$PYTHONPATH:/vsc-hard-mounts/leuven-data/350/vsc35085/programs/genomics_general

#REFNAME=dudPrim
REFNAME=nieuPrim

chrom=1

echo "================="


#9
pop1=(B_NIE_T F_GUE_T P_AVE_T S_HUE_T B_HEI_2000 E_SEV_T)
pop2=(B_DUD_S F_GUE_S P_AVE_S S_COT_S B_HEI_2018 F_CAM_S)


#python /vsc-hard-mounts/leuven-data/350/vsc35085/programs/genomics_general/popgenWindows.py -w 50000 -s 50000 -m 500 --minData 0.1 \
#-g /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$chrom.calls.gz \
#-o /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$chrom.stats_$(echo "${pop1[ID]}")_$(echo "${pop2[ID]}")_w50000_s50000.csv.gz \
#-f phased -T 8 \
#-p $(echo "${pop1[ID]}") \
#-p $(echo "${pop2[ID]}") \
#--popsFile /vsc-hard-mounts/leuven-data/350/vsc35085/scripts/Pogonus_pops_batchALL.txt 

python /vsc-hard-mounts/leuven-data/350/vsc35085/programs/genomics_general/popgenWindows_egglib.py -w 50000 -s 50000 --minSites 1000 --maxMissing 0.25 \
-T 10 --windType coordinate -f phased \
-g /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$chrom.H.calls.gz \
--popsFile /vsc-hard-mounts/leuven-data/350/vsc35085/scripts/Pogonus_pops_batchALL.txt \
-o /lustre1/scratch/350/vsc35085/Maria/Pogonus_reseqALL_$REFNAME.chr_$chrom.stats_$(echo "${pop1[ID]}")_$(echo "${pop2[ID]}")_w50000_s50000_eggStats.stats \
-p $(echo "${pop1[ID]}") \
-p $(echo "${pop2[ID]}") \
-eggB FstWC,Dxy -eggW S,Pi,thetaW,D

#Script written by Steven Van Belleghem
