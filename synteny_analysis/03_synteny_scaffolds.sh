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


#synteny analyses
#SCAFFOLDS A 
ntSynt scaffold_A_dud.fa scaffold_A_nieu.fa -p scaffold_A -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv scaffold_A.synteny_blocks.tsv --fai scaffold_A_dud.fa.fai scaffold_A_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks scaffold_A.synteny_blocks.tsv --sort_order  scaffold_A_dud.fa.fai scaffold_A_nieu.fa.fai --fais > scaffold_A.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai scaffold_A_dud.fa.fai scaffold_A_nieu.fa.fai --prefix scaffold_A --blocks scaffold_A.synteny_blocks.sorted.tsv --length 100 --colour scaffold_A_dud.fa
cp scaffold_A.links.tsv $VSC_DATA
cp scaffold_A.sequence_lengths.tsv $VSC_DATA

#SCAFFOLDS B
ntSynt scaffold_B_dud.fa scaffold_B_nieu.fa -p scaffold_B -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv scaffold_B.synteny_blocks.tsv --fai scaffold_B_dud.fa.fai scaffold_B_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks scaffold_B.synteny_blocks.tsv --sort_order  scaffold_B_dud.fa.fai scaffold_B_nieu.fa.fai --fais > scaffold_B.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai scaffold_B_dud.fa.fai scaffold_B_nieu.fa.fai --prefix scaffold_B --blocks scaffold_B.synteny_blocks.sorted.tsv --length 100 --colour scaffold_B_dud.fa
cp scaffold_B.links.tsv $VSC_DATA
cp scaffold_B.sequence_lengths.tsv $VSC_DATA

#SCAFFOLDS C
ntSynt scaffold_C_dud.fa scaffold_C_nieu.fa -p scaffold_C -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv scaffold_C.synteny_blocks.tsv --fai scaffold_C_dud.fa.fai scaffold_C_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks scaffold_C.synteny_blocks.tsv --sort_order  scaffold_C_dud.fa.fai scaffold_C_nieu.fa.fai --fais > scaffold_C.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai scaffold_C_dud.fa.fai scaffold_C_nieu.fa.fai --prefix scaffold_C --blocks scaffold_C.synteny_blocks.sorted.tsv --length 100 --colour scaffold_C_dud.fa
cp scaffold_C.links.tsv $VSC_DATA
cp scaffold_C.sequence_lengths.tsv $VSC_DATA

#SCAFFOLDS D
ntSynt scaffold_D_dud.fa scaffold_D_nieu.fa -p scaffold_D -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv scaffold_D.synteny_blocks.tsv --fai scaffold_D_dud.fa.fai scaffold_D_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks scaffold_D.synteny_blocks.tsv --sort_order  scaffold_D_dud.fa.fai scaffold_D_nieu.fa.fai --fais > scaffold_D.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai scaffold_D_dud.fa.fai scaffold_D_nieu.fa.fai --prefix scaffold_D --blocks scaffold_D.synteny_blocks.sorted.tsv --length 100 --colour scaffold_D_dud.fa
cp scaffold_D.links.tsv $VSC_DATA
cp scaffold_D.sequence_lengths.tsv $VSC_DATA

#SCAFFOLDS E
ntSynt scaffold_E_dud.fa scaffold_E_nieu.fa -p scaffold_E -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv scaffold_E.synteny_blocks.tsv --fai scaffold_E_dud.fa.fai scaffold_E_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks scaffold_E.synteny_blocks.tsv --sort_order  scaffold_E_dud.fa.fai scaffold_E_nieu.fa.fai --fais > scaffold_E.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai scaffold_E_dud.fa.fai scaffold_E_nieu.fa.fai --prefix scaffold_E --blocks scaffold_E.synteny_blocks.sorted.tsv --length 100 --colour scaffold_E_dud.fa
cp scaffold_E.links.tsv $VSC_DATA
cp scaffold_E.sequence_lengths.tsv $VSC_DATA

#SCAFFOLDS F
ntSynt scaffold_F_dud.fa scaffold_F_nieu.fa -p scaffold_F -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv scaffold_F.synteny_blocks.tsv --fai scaffold_F_dud.fa.fai scaffold_F_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks scaffold_F.synteny_blocks.tsv --sort_order  scaffold_F_dud.fa.fai scaffold_F_nieu.fa.fai --fais > scaffold_F.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai scaffold_F_dud.fa.fai scaffold_F_nieu.fa.fai --prefix scaffold_F --blocks scaffold_F.synteny_blocks.sorted.tsv --length 100 --colour scaffold_F_dud.fa
cp scaffold_F.links.tsv $VSC_DATA
cp scaffold_F.sequence_lengths.tsv $VSC_DATA

#SCAFFOLDS G
ntSynt scaffold_G_dud.fa scaffold_G_nieu.fa -p scaffold_G -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv scaffold_G.synteny_blocks.tsv --fai scaffold_G_dud.fa.fai scaffold_G_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks scaffold_G.synteny_blocks.tsv --sort_order  scaffold_G_dud.fa.fai scaffold_G_nieu.fa.fai --fais > scaffold_G.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai scaffold_G_dud.fa.fai scaffold_G_nieu.fa.fai --prefix scaffold_G --blocks scaffold_G.synteny_blocks.sorted.tsv --length 100 --colour scaffold_G_dud.fa
cp scaffold_G.links.tsv $VSC_DATA
cp scaffold_G.sequence_lengths.tsv $VSC_DATA

#SCAFFOLDS H
ntSynt scaffold_H_dud.fa scaffold_H_nieu.fa -p scaffold_H -t 24 -d 0.01
python denovo_synteny_block_stats.py --tsv scaffold_H.synteny_blocks.tsv --fai scaffold_H_dud.fa.fai scaffold_H_nieu.fa.fai
python sort_ntsynt_blocks.py --synteny_blocks scaffold_H.synteny_blocks.tsv --sort_order  scaffold_H_dud.fa.fai scaffold_H_nieu.fa.fai --fais > scaffold_H.synteny_blocks.sorted.tsv
python format_blocks_gggenomes.py --fai scaffold_H_dud.fa.fai scaffold_H_nieu.fa.fai --prefix scaffold_H --blocks scaffold_H.synteny_blocks.sorted.tsv --length 100 --colour scaffold_H_dud.fa
cp scaffold_H.links.tsv $VSC_DATA
cp scaffold_H.sequence_lengths.tsv $VSC_DATA
