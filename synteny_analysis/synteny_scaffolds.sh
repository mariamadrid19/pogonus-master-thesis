

#SCAFFOLDS A 
awk 'NR==2 {print substr($0, 9100059, 39608137-9100059+1)}' LG01_dud.fasta > scaffold_A_dud.fa #Dudzele 2
sed -i '1i >scaffold_2_RagTag' scaffold_A_dud.fa

awk 'NR==2 {print substr($0, 54672198, 93883070-54672198+1)}' LG01_nieu.fasta > scaffold_A_nieu.fa #Nieuwpoort 3
sed -i '1i >scaffold_3_RagTag' scaffold_A_nieu.fa

#SCAFFOLDS B
awk 'NR==2 {print substr($0, 74064473, 103390037-74064473+1)}' LG01_dud.fasta > scaffold_B_dud.fa
sed -i '1i >scaffold_4_RagTag' scaffold_B_dud.fa

awk 'NR==2 {print substr($0, 120572473, 153501276-120572473+1)}' LG01_nieu.fasta > scaffold_B_nieu.fa
sed -i '1i >scaffold_7_RagTag' ext_scaffold_B_prim_nieu.fa
