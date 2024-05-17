#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name synteny_scaffolds
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -o synteny_scaffolds.%j.out
#SBATCH -A lp_svbelleghem

#SCAFFOLDS A 
awk 'NR==2 {print substr($0, 9100059, 39608137-9100059+1)}' LG01_dud.fasta > scaffold_A_dud.fa #Dudzele 2
sed -i '1i >scaffold_2_RagTag' scaffold_A_dud.fa

awk 'NR==2 {print substr($0, 54672198, 93883070-54672198+1)}' LG01_nieu.fasta > scaffold_A_nieu.fa #Nieuwpoort 3
sed -i '1i >scaffold_3_RagTag' scaffold_A_nieu.fa

#SCAFFOLDS B
awk 'NR==2 {print substr($0, 74064473, 103390037-74064473+1)}' LG01_dud.fasta > scaffold_B_dud.fa #Dudzele 4
sed -i '1i >scaffold_4_RagTag' scaffold_B_dud.fa

awk 'NR==2 {print substr($0, 120572473, 153501276-120572473+1)}' LG01_nieu.fasta > scaffold_B_nieu.fa #Nieuwpoort 7
sed -i '1i >scaffold_7_RagTag' ext_scaffold_B_prim_nieu.fa

#SCAFFOLDS C
awk 'NR==2 {print substr($0, 3681807, 34037187-3681807+1)}' LG02_dud.fasta > scaffold_C_dud.fa #Dudzele 3
sed -i '1i >scaffold_3_RagTag' scaffold_C_dud.fa

awk 'NR==2 {print substr($0, 7103728, 45349787-7103728+1)}' LG02_nieu.fasta > scaffold_C_nieu.fa #Nieuwpoort 4
sed -i '1i >scaffold_4_RagTag' scaffold_C_nieu.fa

#SCAFFOLDS D
awk 'NR==2 {print substr($0, 34117388, 72878678-34117388+1)}' LG02_dud.fasta > scaffold_D_dud.fa #Dudzele 1
sed -i '1i >scaffold_1_RagTag' scaffold_D_dud.fa

awk 'NR==2 {print substr($0, 45472388, 79512942-45472388+1)}' LG02_nieu.fasta > scaffold_D_nieu.fa #Nieuwpoort 5
sed -i '1i >scaffold_5_RagTag' scaffold_D_nieu.fa

#SCAFFOLDS E
awk 'NR==2 {print substr($0, 22340581, 46000200-22340581+1)}' LG03_dud.fasta > scaffold_E_dud.fa #Dudzele 7
sed -i '1i >scaffold_7_RagTag' scaffold_E_dud.fa

awk 'NR==2 {print substr($0, 43076330, 84460820-43076330+1)}' LG03_nieu.fasta > scaffold_E_nieu.fa #Nieuwpoort 2
sed -i '1i >scaffold_2_RagTag' scaffold_E_nieu.fa

#SCAFFOLDS F 
awk 'NR==2 {print substr($0, 13335419, 33377731-13335419+1)}' LG04_dud.fasta > scaffold_F_dud.fa #Dudzele 9
sed -i '1i >scaffold_9_RagTag' scaffold_F_dud.fa

awk 'NR==2 {print substr($0, 24298258, 52146979-24298258+1)}' LG04_nieu.fasta > scaffold_F_nieu.fa  #Nieuwpoort 12
sed -i '1i >scaffold_12_RagTag' scaffold_F_nieu.fa

#SCAFFOLDS G
awk 'NR==2 {print substr($0, 24993515, 44107433-24993515+1)}' LG05_dud.fasta > scaffold_G_dud.fa #Dudzele 11
sed -i '1i >scaffold_11_RagTag' scaffold_G_dud.fa

awk 'NR==2 {print substr($0, 26270178, 54329264-26270178+1)}' LG05_nieu.fasta > scaffold_G_nieu.fa #Nieuwpoort 10
sed -i '1i >scaffold_10_RagTag' scaffold_G_nieu.fa

#SCAFFOLDS H
awk 'NR==2 {print substr($0, 3638244, 23453803-3638244+1)}' LG10_dud.fasta > scaffold_H_dud.fa #Dudzele 10
sed -i '1i >scaffold_10_RagTag' scaffold_H_dud.fa

awk 'NR==2 {print substr($0, 957591, 23831487-957591+1)}' LG10_nieu.fasta > scaffold_H_nieu.fa #Nieuwpoort 17
sed -i '1i >scaffold_17_RagTag' scaffold_H_nieu.fa




