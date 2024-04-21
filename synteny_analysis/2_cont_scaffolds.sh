#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name cont_scaffolds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -o cont_scaffolds.%j.out
#SBATCH -A lp_svbelleghem

#Now it's necessary to study the synteny between the most contiguous (i.e. largest and longest) scaffolds. This is done by observing the AGP files from RagTag and extracting the largest
#scaffolds that match with the Fst signals from the popgen analysis

#To extract the linkage groups from the primary assemblies
awk '/^>LG01_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG01_dud.fasta
awk '/^>LG02_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG02_dud.fasta
awk '/^>LG03_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG03_dud.fasta
awk '/^>LG04_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG04_dud.fasta
awk '/^>LG05_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG05_dud.fasta
awk '/^>LG06_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG06_dud.fasta
awk '/^>LG07_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG07_dud.fasta
awk '/^>LG08_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG08_dud.fasta
awk '/^>LG09_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG09_dud.fasta
awk '/^>LG10_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_dud.fasta > LG10_dud.fasta

awk '/^>LG01_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG01_nieu.fasta
awk '/^>LG02_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG02_nieu.fasta
awk '/^>LG03_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG03_nieu.fasta
awk '/^>LG04_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG04_nieu.fasta
awk '/^>LG05_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG05_nieu.fasta
awk '/^>LG06_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG06_nieu.fasta
awk '/^>LG07_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG07_nieu.fasta
awk '/^>LG08_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG08_nieu.fasta
awk '/^>LG09_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG09_nieu.fasta
awk '/^>LG10_P_chalceus_RagTag$/{flag=1; print; next} flag && /^>/{flag=0} flag && NR>1' 10_sorted_prim_nieu.fasta > LG10_nieu.fasta

#The to extract specific scaffolds (given the positions found in the corresponding AGP files). It's necessary to add a new heading since awk only extracts the raw sequence

#SCAFFOLD 1

#SCAFFOLD 2
awk 'NR==2 {print substr($0, 9100059, 39608137-9100059+1)}' LG01_dud.fasta > ext_scaffold_2_prim_dud.fa
sed -i '1i >scaffold_2_RagTag' ext_scaffold_2_prim_dud.fa

awk 'NR==2 {print substr($0, 43076330, 84460820-43076330+1)}' LG03_dud.fasta > ext_scaffold_2_prim_nieu.fa
sed -i '1i >scaffold_2_RagTag' ext_scaffold_2_prim_nieu.fa

#SCAFFOLD 3

#SCAFFOLD 4

#SCAFFOLD 7

#SCAFFOLD 8

#SCAFFOLD 9

#SCAFFOLD 11
